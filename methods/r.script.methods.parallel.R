################################################################################
##  Name:         r.script.methods.parallel.R
##  Created:      2018.08.02
##  Last edited:  2018.08.02
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for setting up the workers
################################################################################
fn.set.up.worker <- function(ls.cluster, Estim.Opt){
    ##  Load packages to worker
    clusterCall(ls.cluster, function(str.p){
        lapply(str.p, require, character.only = TRUE)
        NULL ## Avoid sending unneccesary info to master
    }, Estim.Opt$str.packages)
    
    ##  Set seed on worker
    clusterSetRNGStream(cl = ls.cluster, iseed = Estim.Opt$i.seed)
    
    ##  Set up functions on the worker and load Estim.Opt
    clusterExport(ls.cluster, c("fn.make.beta.names",
                                "fn.make.beta",
                                "fn.log.lik",
                                "ls.str.par.names",
                                "Estim.Opt"))
    
    ##  Export the row index to the worker
    clusterApply(ls.cluster, ls.data.index, function(ls.tmp){
        ls.index <<- ls.tmp
        NULL
    })
    rm(ls.data.index, envir = .GlobalEnv)
    
    ##  Export the attribute data to the worker
    clusterApply(ls.cluster, ls.data.X, function(ls.tmp){
        ls.X <<- ls.tmp
        NULL
    })
    rm(ls.data.X, envir = .GlobalEnv)
    
    ##  Export the choice data to the worker
    clusterApply(ls.cluster, ls.data.Y, function(ls.tmp){
        ls.Y <<- ls.tmp
        NULL
    })
    rm(ls.data.Y, envir = .GlobalEnv)
    
    ##  Export the interactions data
    if(length(Estim.Opt$ls.het.par) > 0L){
        clusterApply(ls.cluster, ls.data.H, function(m.tmp){
            m.H <<- m.tmp
            NULL
        })
        rm(ls.data.H, envir = .GlobalEnv)
    }
    
    ##  Export the relative scale data
    if(Estim.Opt$b.relative.scale){
        clusterApply(ls.cluster, ls.data.R, function(m.tmp){
            m.R <<- m.tmp
            NULL
        })
        rm(ls.data.R, envir = .GlobalEnv)
    }
    
    ##  Export the class probability data
    if(Estim.Opt$b.latent.class){
        clusterApply(ls.cluster, ls.data.C, function(ls.tmp){
            ls.C <<- ls.tmp
            NULL
        })
        rm(ls.data.C, envir = .GlobalEnv)
        
        ##  Export the list/matrix of constraints
        if(Estim.Opt$b.equality.constrained){
            if(length(Estim.Opt$ls.rand.par) > 0L){
                clusterExport(ls.cluster, "ls.constraints")
                rm(ls.constraints, envir = .GlobalEnv)
            } else {
                clusterExport(ls.cluster, "m.constraints")
                rm(m.constraints, envir = .GlobalEnv)
            }
        }
    }
    
    ##  Export the random draws
    if(Estim.Opt$b.make.draws){
        clusterApply(ls.cluster, ls.draws, function(m.tmp){
            m.draws <<- m.tmp
            NULL
        })
        rm(ls.draws, envir = .GlobalEnv)
        
        clusterApply(ls.cluster, ls.draws.index, function(ls.tmp){
            ls.draws.index <<- ls.tmp
            NULL
        })
        rm(ls.draws.index, envir = .GlobalEnv)
    }
    
    ##  Export a "core counter"
    clusterApply(ls.cluster, as.list(1L:Estim.Opt$i.cores), function(i.x){
        i.am.core.number <<- i.x
        NULL
    })
    
}

################################################################################
##  Function for checking the worker
################################################################################
fn.check.worker <- function(ls.cluster, Estim.Opt){
    ##  Get the objects in .GlobalEnv on the worker
    ls.GlobalEnv.worker <- parLapply(ls.cluster, seq_along(ls.cluster),
                                     function(cl.x){
                                         ls.tmp <- ls(.GlobalEnv)
                                         return(ls.tmp)
                                     })
    
    ##  Get the size of the objects in .GlobalEnv on worker
    ls.ObjectSize.worker <- parLapply(ls.cluster, seq_along(ls.cluster),
                                      function(cl.x){
                                          ls.tmp <- lapply(ls(.GlobalEnv), function(str.x){
                                             i.tmp <- object_size(get(str.x, envir = .GlobalEnv))
                                             return(i.tmp)
                                          })
                                      })
    
    ##  Combine the lists to a matrix
    ls.objects.worker <- mapply(cbind, ls.GlobalEnv.worker, ls.ObjectSize.worker,
                                SIMPLIFY = FALSE)
    
    ls.objects.worker <- lapply(ls.objects.worker, function(m.x){
        colnames(m.x) <- c("object", "size (bytes)")
        return(m.x)
    })
    
    ##  Get a list of the packages loaded on the worker
    ls.packages.worker <- parLapply(ls.cluster, seq_along(ls.cluster),
                                    function(cl.x){
                                        return(search())
                                    })
    
    ##  Set up a list with length equal to the number of workers
    ls.workers <- vector(mode = "list", length = Estim.Opt$i.cores)
    names(ls.workers) <- paste0("worker.", seq_len(Estim.Opt$i.cores))
    
    for(i.i in seq_len(Estim.Opt$i.cores)){
        ls.workers[[i.i]] <- list(ls.objects.worker[[i.i]],
                                  ls.packages.worker[[i.i]])
    }
    ls.workers <- lapply(ls.workers, function(ls.x){
        names(ls.x) <- c("objects", "packages")
        return(ls.x)
    })
    
    ##  Print information about the worker
    cat("###################################################################\n")
    cat("Printing information on workers.")
    print(ls.workers)
    cat("###################################################################\n")
    
    
    ############################################################################
    rm(ls.GlobalEnv.worker, ls.ObjectSize.worker, ls.objects.worker,
       ls.packages.worker, ls.workers)
    ############################################################################
}
