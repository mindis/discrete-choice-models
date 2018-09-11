################################################################################
##  Name:         r.script.methods.parallel.R
##  Created:      2018.08.02
##  Last edited:  2018.09.10
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for setting up the workers
################################################################################
fnSetUpWorker <- function(lsCluster, EstimOpt){
    ##  Load packages to worker
    clusterCall(lsCluster, function(strPack){
        lapply(strPack, require, character.only = TRUE)
        NULL ## Avoid sending unneccesary info to master
    }, EstimOpt$strPackages)
    
    ##  Set seed on worker
    clusterSetRNGStream(cl = lsCluster, iseed = EstimOpt$iSeed)
    
    ##  Set up functions on the worker and load EstimOpt
    clusterExport(lsCluster, c("fnMakeBetaNames",
                               "fnMakeBeta",
                               "fnLogLik",
                               "lsParNames",
                               "EstimOpt"))
    
    ##  Export the row index to the worker
    # clusterApply(lsCluster, ls.data.index, function(ls_tmp){
    #     ls.index <<- ls_tmp
    #     NULL
    # })
    # rm(ls.data.index, envir = .GlobalEnv)
    
    ##  Export the attribute data to the worker
    clusterApply(lsCluster, lsDataX, function(ls_tmp){
        lsX <<- ls_tmp
        NULL
    })
    rm(lsDataX, envir = .GlobalEnv)
    
    ##  Export the choice data to the worker
    clusterApply(lsCluster, lsDataY, function(ls_tmp){
        lsY <<- ls_tmp
        NULL
    })
    rm(lsDataY, envir = .GlobalEnv)
    
    ##  Export the interactions data
    if(length(EstimOpt$lsP_het) > 0L){
        clusterApply(lsCluster, lsDataH, function(m_tmp){
            mH <<- m_tmp
            NULL
        })
        rm(lsDataH, envir = .GlobalEnv)
    }
    
    ##  Export the relative scale data
    if(EstimOpt$bRelativeScale){
        clusterApply(lsCluster, lsDataR, function(m_tmp){
            mR <<- m_tmp
            NULL
        })
        rm(lsDataR, envir = .GlobalEnv)
    }
    
    ##  Export the class probability data
    if(EstimOpt$bLatentClass){
        clusterApply(lsCluster, lsDataC, function(ls_tmp){
            lsC <<- ls_tmp
            NULL
        })
        rm(lsDataC, envir = .GlobalEnv)
        
        ##  Export the list/matrix of constraints
        if(EstimOpt$bEqualityConstrained){
            clusterExport(lsCluster, "lsDdelta")
            clusterExport(lsCluster, "mDeltaExpanded")
            rm(lsDelta, mDeltExpanded, envir = .GlobalEnv)
        }
    }
    
    ##  Export the random draws
    if(EstimOpt$bMakeDraws){
        clusterApply(lsCluster, lsDraws, function(m_tmp){
            mDraws <<- m_tmp
            NULL
        })
        rm(lsDraws, envir = .GlobalEnv)
        
        # clusterApply(lsCluster, ls.draws.index, function(ls_tmp){
        #     ls.draws.index <<- ls_tmp
        #     NULL
        # })
        # rm(ls.draws.index, envir = .GlobalEnv)
    }
    
    ##  Export a "core counter"
    clusterApply(lsCluster, as.list(1L:EstimOpt$iCores), function(iX){
        iAmCoreNumber <<- iX
        NULL
    })
    
}

################################################################################
##  Function for checking the worker
################################################################################
fnCheckWorker <- function(lsCluster, EstimOpt){
    ##  Get the objects in .GlobalEnv on the worker
    lsGlobalEnv_worker <- parLapply(lsCluster, seq_along(lsCluster),
                                    function(clX){
                                        ls_tmp <- ls(.GlobalEnv)
                                        return(ls_tmp)
                                    })
    
    ##  Get the size of the objects in .GlobalEnv on worker
    lsObjectSize_worker <- parLapply(lsCluster, seq_along(lsCluster),
                                     function(clX){
                                         ls_tmp <- lapply(ls(.GlobalEnv), function(strX){
                                             dSize <- object_size(get(strX, envir = .GlobalEnv))
                                             return(dSize)
                                         })
                                     })
    
    ##  Combine the lists to a matrix
    lsObjects_worker <- mapply(cbind, lsGlobalEnv_worker, lsObjectSize_worker,
                               SIMPLIFY = FALSE)
    
    lsObjects_worker <- lapply(lsObjects_worker, function(mX){
        colnames(mX) <- c("Object", "Size (bytes)")
        return(mX)
    })
    
    ##  Get a list of the packages loaded on the worker
    lsPackages_worker <- parLapply(lsCluster, seq_along(lsCluster),
                                   function(clX){
                                       return(search())
                                   })
    
    ##  Set up a list with length equal to the number of workers
    lsWorkers <- vector(mode = "list", length = EstimOpt$iCores)
    names(lsWorkers) <- paste0("worker.", seq_len(EstimOpt$iCores))
    
    for(i in seq_len(EstimOpt$iCores)){
        lsWorkers[[i]] <- list(lsObjects_worker[[i]],
                               lsPackages_worker[[i]])
    }
    lsWorkers <- lapply(lsWorkers, function(lsX){
        names(lsX) <- c("Objects", "Packages")
        return(lsX)
    })
    
    ##  Print information about the worker
    cat("###################################################################\n")
    cat("Printing information on workers.")
    print(lsWorkers)
    cat("###################################################################\n")
    
    
    ############################################################################
    rm(lsGlobalEnv_worker, lsObjectSize_worker, lsObjects_worker,
       lsPackages_worker, lsWorkers)
    ############################################################################
}