################################################################################
##  Name:         r.script.methods.run.R
##  Created:      2018.08.02
##  Last edited:  2018.08.02
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for searching for starting values
################################################################################
fn.search.starting.values <- function(v.param, v.str.par.names, ls.cluster,
                                      Estim.Opt){
    ##  Check that we have specified sufficient number of models and correct
    if(Estim.Opt$i.nr.of.models > Estim.Opt$i.nr.of.starting.models){
        Estim.Opt$i.nr.of.starting.models <<- Estim.Opt$i.nr.of.models * 2
    }

    ##  Set up a matrix of v.param repeated equal to nr of starting models
    m.param <- matrix(rep(v.param, Estim.Opt$i.nr.of.starting.models),
                      ncol = length(v.param), byrow = TRUE) + 0.001
    
    ##  Set up a matrix of random numbers to use in calculation
    m.runif <- matrix(runif(Estim.Opt$i.nr.of.starting.models * length(v.param)),
                      ncol = length(v.param)) - 0.5
    
    ##  Create the matrix of starting values
    m.param <- m.param + (m.param * m.runif * Estim.Opt$d.multiplier * 2)
    
    ##  Turn into a list for faster processing
    ls.param <- as.list(data.frame(t(m.param)))
    ls.param <- lapply(ls.param, function(v.x){
        names(v.x) <- v.str.par.names
        return(v.x)
    })
    
    ##  Check whether we are estimating models in parallel
    if(Estim.Opt$b.parallel){
        ls.ll.values <- lapply(ls.param, function(v.x){
            ls.ll <- parLapply(ls.cluster, seq_along(ls.cluster),
                               function(i.am.not.using.this){
                                   fn.log.lik(v.x)
                               })
            return(sum(Reduce(c, ls.ll)))
        })
    } else{
        ls.ll.values <- lapply(ls.param, function(v.x){
            sum(fn.log.lik(v.x))
        })
    }
    
    ##  Create a matrix and sort based on LL values
    m.param <- cbind(m.param, Reduce(rbind, ls.ll.values))
    m.param <- m.param[order(m.param[, ncol(m.param)], decreasing = TRUE), ]
    return(m.param[1L:Estim.Opt$i.nr.of.models, 1L:length(v.param), drop = FALSE])
}

################################################################################
##  Function for running the model
################################################################################
fn.run.model <- function(Estim.Opt){
    ##  Check that file writing is turned oof
    if(sink.number() > 0L) sink()
    
    ##  Set seed
    set.seed(Estim.Opt$i.seed)
    
    ##  If we are not using parallel computing, set number of cores to 1
    if(!Estim.Opt$b.parallel) Estim.Opt$i.cores <- 1L
    
    ##  Check number of cores
    if(Estim.Opt$b.parallel && Estim.Opt$i.cores >= detectCores()){
        Estim.Opt$i.cores <- max(1L, detectCores() - 1L)
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
    fn.set.up.data(Estim.Opt)
    time.tmp <- fn.time(proc.time() - proc.time.start)
    cat("Setting up the data took",
        time.tmp[[1]], "days",
        time.tmp[[2]],"hours",
        time.tmp[[3]], "minutes and",
        time.tmp[[4]], "seconds.\n")
    cat("###################################################################\n")
    cat("\n")
    rm(time.tmp)
    
    ##  Generate parameter names
    ls.str.par.names <<- fn.make.beta.names(Estim.Opt)
    v.str.par.names <- Reduce(c, ls.str.par.names)
    
    ############################################################################
    ##  Check whether we are estimating on a single or multiple cores
    ############################################################################
    if(Estim.Opt$b.parallel){
        ##  Set up the cluster
        ls.cluster <- makeCluster(Estim.Opt$i.cores, type = "PSOCK",
                                  manual = FALSE,
                                  outfile = Estim.Opt$str.output.debug)
        
        ##  Make sure that the cluster is closed when the function exits
        on.exit(stopCluster(ls.cluster))
        
        ##  Set up the workers -- cleaning is done within this function
        fn.set.up.worker(ls.cluster, Estim.Opt)
        
        ##  Check the worker -- functionality to be added
        if(Estim.Opt$b.print.worker.info){
            fn.check.worker(ls.cluster, Estim.Opt)
        }
        
        ##  Garbage collection
        invisible(gc(verbose = FALSE))
    } else {
        ##  Save the data witih the correct names globally
        ls.index <<- ls.data.index[[1L]]
        rm(ls.data.index, envir = .GlobalEnv)
        
        ls.X <<- ls.data.X[[1L]]
        rm(ls.data.X, envir = .GlobalEnv)
        
        ls.Y <<- ls.data.Y[[1L]]
        rm(ls.data.Y, envir = .GlobalEnv)
        
        if(length(Estim.Opt$ls.het.par) > 0L){
            m.H <<- ls.data.H[[1L]]
            rm(ls.data.H, envir = .GlobalEnv)
        } 
        
        if(Estim.Opt$b.relative.scale){
            m.R <<- ls.data.R[[1L]]
            rm(ls.data.R, envir = .GlobalEnv)
        } 
        
        if(Estim.Opt$b.latent.class){
            ls.C <<- ls.data.C[[1L]]
            rm(ls.data.C, envir = .GlobalEnv)
        } 
        
        if(Estim.Opt$b.make.draws){
            m.draws <<- ls.draws[[1L]]
            rm(ls.draws, envir = .GlobalEnv)
            ls.draws.index <<- ls.draws.index[[1L]]
        }
        
        ##  Set empty ls.cluster for starting value search
        ls.cluster <- NULL
        
        ##  Garbage collection
        invisible(gc(verbose = FALSE))
    }
    
    ############################################################################
    ##  Check whether we are searching for starting values
    ############################################################################
    if(Estim.Opt$b.search.starting.values && Estim.Opt$i.iterlim > 0L){
        proc.time.start <- proc.time()
        cat("###################################################################\n")
        cat("Searching", Estim.Opt$i.nr.of.starting.models,
            "models for starting values ... \n")
        m.param <- fn.search.starting.values(v.param, v.str.par.names, ls.cluster,
                                             Estim.Opt)
        time.tmp <- fn.time(proc.time() - proc.time.start)
        cat("Starting value search took",
            time.tmp[[1]], "days",
            time.tmp[[2]],"hours",
            time.tmp[[3]], "minutes and",
            time.tmp[[4]], "seconds.\n")
        cat("Now estimating", Estim.Opt$i.nr.of.models, "model(s) to completion!\n")
        cat("###################################################################\n")
        cat("\n")
        rm(time.tmp)
    }
    
    ############################################################################
    ##  Write the wrapper function for the log-likelihood
    ############################################################################
    fn.log.lik.wrapper <- function(v.param){
        ##  Check whether we are using serial or parallel processing
        if(Estim.Opt$b.parallel){
            ls.log.lik <- parLapply(ls.cluster, seq_along(ls.cluster),
                                    function(i.am.not.using.this){
                                        fn.log.lik(v.param)
                                    })
            v.log.lik <- Reduce(c, ls.log.lik)
        } else {
            v.log.lik <- fn.log.lik(v.param)
        }
        
        ##  Print progress bar for hessian matrix calculation
        if(exists("ls.model.tmp", envir = environment())){
            i.count <<- i.count + 1L
            i.K <- length(v.param)
            i.call.numDeriv <- 4L * i.K^2L + 4L * i.K + 3L
            v.percentile <- floor(seq(0L, i.call.numDeriv, by = (i.call.numDeriv / 10)))
            if(i.count == 1L){cat("Hessian calculation: 0% ...")}
            if(i.count == v.percentile[2L]){cat("10% ...")}
            if(i.count == v.percentile[3L]){cat("20% ...")}
            if(i.count == v.percentile[4L]){cat("30% ...")}
            if(i.count == v.percentile[5L]){cat("40% ...")}
            if(i.count == v.percentile[6L]){cat("50% ...")}
            if(i.count == v.percentile[7L]){cat("60% ...")}
            if(i.count == v.percentile[8L]){cat("70% ...")}
            if(i.count == v.percentile[9L]){cat("80% ...")}
            if(i.count == v.percentile[10L]){cat("90% ...")}
            if(i.count == v.percentile[11L]){cat("100% ...")}
        }
        return(v.log.lik)
    }
    
    ############################################################################
    ##  Loop over the number of models to run to completion
    ############################################################################
    i.count.failed <<- 0L
    if(!Estim.Opt$b.search.starting.values) Estim.Opt$i.nr.of.models <- 1L
    for(i.n in seq_len(Estim.Opt$i.nr.of.models)){
        proc.time.start <- proc.time()
        
        cat("\n")
        cat("###################################################################\n")
        cat("The model is estimated using ", Estim.Opt$i.cores,  "core(s).\n")
        cat("Estimating model ", i.n, " out of ", Estim.Opt$i.nr.of.models,
            ". (", i.count.failed, " model(s) have failed).\n", sep = "")
        cat("###################################################################\n")
        cat("\n")
        
        ##  Get the starting values and name them for reference
        if(Estim.Opt$b.search.starting.values && Estim.Opt$i.iterlim > 0){
            v.param <<- m.param[i.n, ]
            names(v.param) <<- v.str.par.names
        } else {
            names(v.param) <<- v.str.par.names
        }
        
        ########################################################################
        ##  Run the model
        ########################################################################
        ls.model.tmp <- maxLik(fn.log.lik.wrapper,
                               start = v.param,
                               method = Estim.Opt$str.method,
                               print.level = Estim.Opt$i.print.level,
                               tol = Estim.Opt$d.tol,
                               reltol = Estim.Opt$d.reltol,
                               gradtol = Estim.Opt$d.gradtol,
                               steptol = Estim.Opt$d.steptol,
                               iterlim = Estim.Opt$i.iterlim,
                               finalHessian = Estim.Opt$b.final.hessian)
        
        ########################################################################
        ##  Check for convergence
        ########################################################################
        if(ls.model.tmp$code == 0L || ls.model.tmp$code == 2L){
            cat("###################################################################\n")
            cat("Calculating the final Hessian -- get some coffee!\n\n")
            ##  Set global counter for the "progress bar"
            i.count <<- 0L
            
            ##  Sum of the wrapper functoin
            fn.log.lik.wrapper.sum <- function(v.param.est){
                sum(fn.log.lik.wrapper(v.param.est))
            }
            
            ##  Calculate the Hessian matrix using numDeriv()
            if(!Estim.Opt$b.final.hessian){
                ls.model.tmp$hessian <- numDeriv::hessian(func = fn.log.lik.wrapper.sum,
                                                          x = ls.model.tmp$estimate)
            }

            cat("\n")
            ##  Approximate the log.lik(0)
            ls.model.tmp$d.ll.zero <- Estim.Opt$i.obs * log(1L/Estim.Opt$i.alts)
            
            ##  Remove the counter from the global environment
            rm(i.count, envir = .GlobalEnv)
            
            ##  Add additional information to the model object
            ls.model.tmp$time <- proc.time() - proc.time.start
            ls.model.tmp$v.starting.values <- v.param
            ls.model.tmp$i.model.number <- i.n
            ls.model.tmp$i.count.failed <- i.count.failed
            ls.model.tmp$i.seed <- Estim.Opt$i.seed
            
            cat("###################################################################\n")
            cat("Estimation of model", i.n, "complete.\n")
            
            ##  Save model object globally and to disk
            assign(paste0("ls.model.", i.n), ls.model.tmp, envir = .GlobalEnv)
            saveRDS(ls.model.tmp,
                    file = paste("MODEL", Estim.Opt$str.output.estimates,
                                 i.n, "rds", sep = "."))
            
            ##  Sink the results to a .txt file for easy inspection
            if(i.n == 1L){
                sink(paste0(Estim.Opt$str.output.estimates, ".txt"))
            } else {
                sink(paste0(Estim.Opt$str.output.estimates, ".txt"),
                     append = TRUE)
            }
            print(fn.model.summary(ls.model.tmp))
            sink()
            
            ##  Print summary to console
            print(fn.model.summary(ls.model.tmp))
        } else {
            cat("###################################################################\n")
            cat("Estimation of model", i.n, "failed.")
            i.count.failed <<- i.count.failed + 1L
        }
        rm(ls.model.tmp)
        invisible(gc(verbose = FALSE))
    }
}
