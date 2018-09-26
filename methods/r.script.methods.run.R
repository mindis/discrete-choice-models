################################################################################
##  Name:         r.script.methods.run.R
##  Created:      2018.08.02
##  Last edited:  2018.09.10
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for searching for starting values
################################################################################
fnSearchStartingValues <- function(vP, strParNames, lsCluster,
                                   EstimOpt){
    ##  Set some overall parameters
    iN <- EstimOpt$iN
    iJ <- EstimOpt$iJ
    iT <- EstimOpt$iT
    iM_search <- EstimOpt$iM_search
    iM_run <- EstimOpt$iM_run
    
    ##  Check that we have specified sufficient number of models and correct
    if(iM_run > iM_search){
        iM_search <<- iM_run * 2
    }
    
    ##  Set up a matrix of vP repeated equal to nr of starting models
    mP <- matrix(rep(vP, iM_search), ncol = length(vP), byrow = T) + 0.001
    
    ##  Set up a matrix of random numbers to use in calculation
    mRunif <- matrix(runif(iM_search * length(vP)), ncol = length(vP)) - 0.5
    
    ##  Create the matrix of starting values
    mP <- mP + (mP * mRunif * EstimOpt$dMulti * 2)
    
    ##  Turn into a list for faster processing
    lsP <- as.list(data.frame(t(mP)))
    lsP <- lapply(lsP, function(vP_tmp){
        names(vP_tmp) <- strParNames
        return(vP_tmp)
    })
    
    ##  Check whether we are estimating models in parallel
    if(EstimOpt$bParallel){
        lsLL_values <- lapply(lsP, function(vP_tmp){
            lsLL <- parLapply(lsCluster, seq_along(lsCluster),
                              function(iClusterCounter){
                                  fnLogLik(vP_tmp)
                              })
            return(sum(Reduce(c, lsLL)))
        })
    } else{
        lsLL_values <- lapply(lsP, function(vP_tmp){
            sum(fnLogLik(vP_tmp))
        })
    }
    
    ##  Create a matrix and sort based on LL values
    mP <- cbind(mP, Reduce(rbind, lsLL_values))
    mP <- mP[order(mP[, ncol(mP)], decreasing = T), ]
    return(mP[1L:iM_run, 1L:length(vP), drop = F])
}

################################################################################
##  Function for running the model
################################################################################
fnRunModel <- function(EstimOpt){
    ##  Check that file writing is turned oof
    if(sink.number() > 0L) sink()
    
    ##  Set some overall parameters
    iN <- EstimOpt$iN
    iJ <- EstimOpt$iJ
    iT <- EstimOpt$iT
    iM_search <- EstimOpt$iM_search
    iM_run <- EstimOpt$iM_run
    
    ##  Set seed
    set.seed(EstimOpt$iSeed)
    
    ##  If we are not using parallel computing, set number of cores to 1
    if(!EstimOpt$bParallel) EstimOpt$iCores <- 1L
    
    ##  Check number of cores
    if(EstimOpt$bParallel && EstimOpt$iCores >= detectCores()){
        EstimOpt$iCores <- max(1L, detectCores() - 1L)
        cat("###################################################################\n")
        cat("The number of cores has been adjusted to detectCores() -1L \n")
        cat("###################################################################\n")
        cat("\n")
    }
    
    ############################################################################
    ##  Set up the data
    ############################################################################
    proc.time.start <- proc.time()
    cat("###################################################################\n")
    cat("Setting up the data and generating random draws.\n")
    fnSetUpData(EstimOpt)
    time_tmp <- fnTime(proc.time() - proc.time.start)
    cat("Setting up the data took",
        time_tmp[[1]], "days",
        time_tmp[[2]],"hours",
        time_tmp[[3]], "minutes and",
        time_tmp[[4]], "seconds.\n")
    cat("###################################################################\n")
    cat("\n")
    ############################################################################
    rm(time_tmp)
    ############################################################################
    
    ##  Generate parameter names
    lsParNames <<- fnMakeBetaNames(EstimOpt)
    strParNames <- Reduce(c, lsParNames)
    
    ############################################################################
    ##  Check whether we are estimating on a single or multiple cores
    ############################################################################
    if(EstimOpt$bParallel){
        ##  Set up the cluster
        lsCluster <- makeCluster(EstimOpt$iCores, type = "PSOCK",
                                 manual = F,
                                 outfile = EstimOpt$strDebugFile)
        
        ##  Make sure that the cluster is closed when the function exits
        on.exit(stopCluster(lsCluster))
        
        ##  Set up the workers -- cleaning is done within this function
        fnSetUpWorker(lsCluster, EstimOpt)
        
        ##  Check the worker -- functionality to be added
        if(EstimOpt$bPrintWorkerInfo){
            fnCheckWorker(lsCluster, EstimOpt)
        }
        
        ########################################################################
        invisible(gc(verbose = F))
        ########################################################################
    } else {
        ##  Save the data witih the correct names globally
        # ls.index <<- ls.data.index[[1L]]
        # rm(ls.data.index, envir = .GlobalEnv)
        
        lsX <<- lsDataX[[1L]]
        rm(lsDataX, envir = .GlobalEnv)
        
        lsY <<- lsDataY[[1L]]
        rm(lsDataY, envir = .GlobalEnv)
        
        if(length(EstimOpt$lsP_het) > 0L){
            mH <<- lsDataH[[1L]]
            rm(lsDataH, envir = .GlobalEnv)
        } 
        
        if(EstimOpt$bRelativeScale){
            mR <<- lsDataR[[1L]]
            rm(lsDataR, envir = .GlobalEnv)
        } 
        
        if(EstimOpt$bLatentClass){
            lsC <<- lsDataC[[1L]]
            rm(lsDataC, envir = .GlobalEnv)
        } 
        
        if(EstimOpt$bMakeDraws){
            mDraws <<- lsDraws[[1L]]
            rm(lsDraws, envir = .GlobalEnv)
            # ls.draws.index <<- ls.draws.index[[1L]]
        }
        
        ##  Set empty lsCluster for starting value search
        lsCluster <- NULL
        
        ##  Garbage collection
        invisible(gc(verbose = F))
    }
    
    ############################################################################
    ##  Check whether we are searching for starting values
    ############################################################################
    if(EstimOpt$bSearchStartingValues && EstimOpt$iIterlim > 0L){
        proc.time.start <- proc.time()
        cat("###################################################################\n")
        cat("Searching", iM_search,
            "models for starting values ... \n")
        mP <- fnSearchStartingValues(vP, strParNames, lsCluster, EstimOpt)
        time_tmp <- fnTime(proc.time() - proc.time.start)
        cat("Starting value search took",
            time_tmp[[1]], "days",
            time_tmp[[2]],"hours",
            time_tmp[[3]], "minutes and",
            time_tmp[[4]], "seconds.\n")
        cat("Now estimating", iM_run, "model(s) to completion!\n")
        cat("###################################################################\n")
        cat("\n")
        ########################################################################
        rm(time_tmp)
        ########################################################################
    }
    
    ############################################################################
    ##  Write the wrapper function for the log-likelihood
    ############################################################################
    fnLogLikWrapper <- function(vP){
        ##  Check whether we are using serial or parallel processing
        if(EstimOpt$bParallel){
            lsLogLik <- parLapply(lsCluster, seq_along(lsCluster),
                                  function(iClusterCounter){
                                      fnLogLik(vP)
                                  })
            vLogLik <- Reduce(c, lsLogLik)
        } else {
            vLogLik <- fnLogLik(vP)
        }
        
        ##  Print progress bar for hessian matrix calculation
        if(exists("lsModel_tmp", envir = environment())){
            iCount <<- iCount + 1L
            iK <- length(vP)
            iCallnumDeriv <- 4L * iK^2L + 4L * iK + 3L
            vPercentile <- floor(seq(0L, iCallnumDeriv, by = (iCallnumDeriv / 10)))
            if(iCount == 1L){cat("Hessian calculation: 0% ...")}
            if(iCount == vPercentile[2L]){cat("10% ...")}
            if(iCount == vPercentile[3L]){cat("20% ...")}
            if(iCount == vPercentile[4L]){cat("30% ...")}
            if(iCount == vPercentile[5L]){cat("40% ...")}
            if(iCount == vPercentile[6L]){cat("50% ...")}
            if(iCount == vPercentile[7L]){cat("60% ...")}
            if(iCount == vPercentile[8L]){cat("70% ...")}
            if(iCount == vPercentile[9L]){cat("80% ...")}
            if(iCount == vPercentile[10L]){cat("90% ...")}
            if(iCount == vPercentile[11L]){cat("100% ...")}
        }
        return(vLogLik)
    }
    
    ############################################################################
    ##  Loop over the number of models to run to completion
    ############################################################################
    iCount_failed <<- 0L
    if(!EstimOpt$bSearchStartingValues) iM_run <- 1L
    for(m in seq_len(iM_run)){
        proc.time.start <- proc.time()
        
        cat("\n")
        cat("###################################################################\n")
        cat("The model is estimated using ", EstimOpt$iCores,  "core(s).\n")
        cat("Estimating model ", m, " out of ", iM_run,
            ". (", iCount_failed, " model(s) have failed).\n", sep = "")
        cat("###################################################################\n")
        cat("\n")
        
        ##  Get the starting values and name them for reference
        if(EstimOpt$bSearchStartingValues && EstimOpt$iIterlim > 0){
            vP <<- mP[m, ]
            names(vP) <<- strParNames
        } else {
            names(vP) <<- strParNames
        }
        
        ########################################################################
        ##  Run the model
        ########################################################################
        lsModel_tmp <- maxLik(fnLogLikWrapper,
                              start = vP,
                              method = EstimOpt$strMethod,
                              print.level = EstimOpt$iPrintLevel,
                              tol = EstimOpt$dTol,
                              reltol = EstimOpt$dReltol,
                              gradtol = EstimOpt$dGradtol,
                              steptol = EstimOpt$dSteptol,
                              iterlim = EstimOpt$iIterlim,
                              finalHessian = EstimOpt$bFinalHessian)
        
        ########################################################################
        ##  Check for convergence
        ########################################################################
        if(lsModel_tmp$code == 0L || lsModel_tmp$code == 2L){
            cat("###################################################################\n")
            cat("Calculating the final Hessian -- get some coffee!\n\n")
            ##  Set global counter for the "progress bar"
            iCount <<- 0L
            
            ##  Sum of the wrapper functoin
            fnLogLikWrapperSum <- function(vP_est){
                sum(fnLogLikWrapper(vP_est))
            }
            
            ##  Calculate the Hessian matrix using numDeriv()
            if(!EstimOpt$bFinalHessian){
                lsModel_tmp$hessian <- numDeriv::hessian(func = fnLogLikWrapperSum,
                                                         x = lsModel_tmp$estimate)
            }
            
            cat("\n")
            ##  Approximate the log.lik(0)
            lsModel_tmp$dLL_zero <- EstimOpt$iN_obs * log(1L/iJ)
            
            ##  Remove the counter from the global environment
            rm(iCount, envir = .GlobalEnv)
            
            ##  Add additional information to the model object
            lsModel_tmp$time <- proc.time() - proc.time.start
            lsModel_tmp$vP_start <- vP
            lsModel_tmp$iM <- m
            lsModel_tmp$iCount_failed <- iCount_failed
            lsModel_tmp$iSeed <- EstimOpt$iSeed
            lsModel_tmp$dConvergenceCriteria <- t(colSums2(lsModel_tmp$gradientObs)) %*% vcov(lsModel_tmp) %*% colSums2(lsModel_tmp$gradientObs)
            
            cat("###################################################################\n")
            cat("Estimation of model", m, "complete.\n")
            
            ##  Save model object globally and to disk
            assign(paste0("lsModel", m), lsModel_tmp, envir = .GlobalEnv)
            saveRDS(lsModel_tmp,
                    file = paste("MODEL", EstimOpt$strOutputFile,
                                 m, "rds", sep = "."))
            
            ##  Sink the results to a .txt file for easy inspection
            if(m == 1L){
                sink(paste0(EstimOpt$strOutputFile, ".txt"))
            } else {
                sink(paste0(EstimOpt$strOutputFile, ".txt"),
                     append = T)
            }
            print(fnModelSummary(lsModel_tmp))
            sink()
            
            ##  Print summary to console
            print(fnModelSummary(lsModel_tmp))
        } else {
            cat("###################################################################\n")
            cat("Estimation of model", m, "failed.")
            iCount_failed <<- iCount_failed + 1L
        }
        ########################################################################
        rm(lsModel_tmp)
        invisible(gc(verbose = F))
        ########################################################################
    }
}