################################################################################
##  Name:         r.script.ll.lc.mixl.R
##  Created:      2018.08.03
##  Last edited:  2018.08.04
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

fn.log.lik <- function(v.param){
    ############################################################################
    ##  Class probability function
    ############################################################################
    ##  Zeros are added for the normalizing class
    v.theta <- c(v.param[ls.str.par.names[["str.class"]]],
                 rep(0L, length(Estim.Opt$str.class.par)))
    
    ls.class.prob <- lapply(seq_along(ls.C), function(i.x){
        i.K <- length(Estim.Opt$str.class.par)
        i.start <- 1L + (i.K * (i.x - 1L))
        i.end <- i.K * i.x
        ##  IND*CT
        exp(crossprod(t(ls.C[[i.x]]), v.theta[i.start:i.end]))
    })
    
    ##  Calculate the probability - Vector IND*CT
    v.class.prob.sum <- Reduce("+", ls.class.prob)
    ls.class.prob <- lapply(ls.class.prob, function(v.x){
        v.prob <- v.x / v.class.prob.sum
        v.prob[is.na(v.prob)] <- 0
        as.vector(v.prob)
    })
    
    ##  IND*CT x CLASSES
    m.class.prob <- Reduce(cbind, ls.class.prob)
    ##  CT x IND*CLASSES
    m.class.prob <- matrix(as.vector(m.class.prob), nrow = Estim.Opt$i.tasks)
    
    ############################################################################
    rm(v.theta, ls.class.prob, v.class.prob.sum)
    ############################################################################
    
    ############################################################################
    ##  Make a list of length equal to i.Q containing the betas IND*DRAWS x NVAR
    ############################################################################
    i.Q <- Estim.Opt$i.classes
    str.tmp <- names(Estim.Opt$ls.rand.par)
    m.beta.r <- matrix(v.param[ls.str.par.names[["str.mean"]]], ncol = i.Q)
    m.sigma.r <- matrix(v.param[ls.str.par.names[["str.std"]]], ncol = i.Q)
    rownames(m.beta.r) <- rownames(m.sigma.r) <- paste0("beta.", str.tmp)
    
    ##  IND*DRAWS x NVAR
    ls.beta <- lapply(seq_len(i.Q), function(i.x){
        i.K <- nrow(m.beta.r)
        i.start <- 1 + (i.x - 1) * i.K
        i.end <- i.x * i.K
        m.beta <- fn.make.beta(Estim.Opt, m.beta.r[, i.x], m.sigma.r[, i.x],
                               m.draws[, i.start:i.end])
        return(m.beta)
    })
    
    ############################################################################
    rm(m.beta.r, m.sigma.r, i.Q, str.tmp)
    ############################################################################
    
    ############################################################################
    ##  Interactions with the means m.H -- IND*DRAWS x NVAR
    ############################################################################
    if(length(Estim.Opt$ls.het.par) > 0){
        if(Estim.Opt$b.class.specific){
            i.tmp <- length(ls.str.par.names[["str.het"]]) / Estim.Opt$i.classes
            ##  Different for each class
            for(j in seq_along(ls.beta)){
                i.start <- 1 + (1 - j) * i.tmp
                i.end <- i.tmp * j
                v.phi <- v.param[ls.str.par.names[["str.het"]][i.start:i.end]]
                for(i in seq_along(Estim.Opt$ls.het.par)){
                    str.tmp <- names(Estim.Opt$ls.het.par[i])
                    v.phi.tmp <- v.phi[grep(str.tmp, names(v.phi))]
                    v.tmp <- crossprod(t(m.H[, Estim.Opt$ls.het.par[[i]]]), v.phi.tmp)
                    if(Estim.Opt$ls.rand.par[[str.tmp]] %in% c("-ln", "ln")){
                        ls.beta[[j]][, str.tmp] <- ls.beta[[j]][, str.tmp] * exp(v.tmp)
                    } else{
                        ls.beta[[j]][, str.tmp] <- ls.beta[[j]][, str.tmp] + v.tmp
                    }
                }
            }
        } else{
            ##  Same for all classes
            v.phi <- v.param[ls.str.par.names[["str.het"]]]
            for(i in seq_along(Estim.Opt$ls.het.par)){
                str.tmp <- names(Estim.Opt$ls.het.par[i])
                ##  This might be problematic with partial matching?
                v.phi.tmp <- v.phi[grep(str.tmp, names(v.phi))]
                ##  Multiply data with the parameters
                v.tmp <- crossprod(t(m.H[, Estim.Opt$ls.het.par[[i]]]), v.phi.tmp)
                ##  Multiply with the betas
                for(j in seq_along(ls.beta)){
                    if(Estim.Opt$ls.rand.par[[str.tmp]] %in% c("-ln", "ln")){
                        ls.beta[[j]][, str.tmp] <- ls.beta[[j]][, str.tmp] * exp(v.tmp)
                    } else{
                        ls.beta[[j]][, str.tmp] <- ls.beta[[j]][, str.tmp] + v.tmp
                    }
                }
            }
        }
    }
    
    ############################################################################
    ##  Check if we are calculating utility in WTP space -- does not work with 
    ##  fixed parameters
    ############################################################################
    if(Estim.Opt$b.wtp.space){
        for(i in seq_along(ls.beta)){
            v.cost <- ls.beta[[i]][, Estim.Opt$str.cost]
            ls.beta[[i]][, Estim.Opt$str.cost] <- 1L
            ##  IND*DRAWS x NVAR
            ls.beta[[i]] <- ls.beta[[i]] * v.cost
        }
    }
    
    ############################################################################
    ##  Calculate the random part of utility
    ############################################################################
    ##  Subset the data to get the variables with random parameters
    ls.X.r <- lapply(ls.X, function(m.x){
        ##  IND*CT x NVAR
        m.x[, names(Estim.Opt$ls.rand.par), drop = FALSE]
    })
    
    i.T <- Estim.Opt$i.tasks
    i.D <- Estim.Opt$i.draws
    i.ind <- nrow(ls.X.r[[1L]]) / i.T
    
    ##  Create list of lists
    ls.utility <- vector(mode = "list", length = Estim.Opt$i.classes)
    ls.utility <- lapply(ls.utility, function(x){
        ls.tmp <- lapply(seq_len(Estim.Opt$i.alts), function(i.x){
            matrix(NA, nrow = nrow(ls.X.r[[1L]]), ncol = Estim.Opt$i.draws)
        })
        return(ls.tmp)
    })
    
    ##  IND*CT x DRAWS
    for(i.c in 1L:Estim.Opt$i.classes){
        for(i.i in 1L:i.ind){
            v.rows <- (1L + ((i.i - 1L) * i.T)):(i.i * i.T)
            for(i.j in 1L:Estim.Opt$i.alts){
                ls.utility[[i.c]][[i.j]][v.rows, ] <- tcrossprod(ls.X.r[[i.j]][v.rows, , drop = FALSE],
                                                          ls.beta[[i.c]][(1L + ((i.i - 1L) * i.D)):(i.i * i.D), ])
            }
        }
    }
    
    ############################################################################
    ##  Calculate the fixed part of utility
    ############################################################################
    if(length(Estim.Opt$str.fixed.par) > 0){
        m.beta.f <- matrix(v.param[ls.str.par.names[["str.fixed"]]],
                           ncol = Estim.Opt$i.classes)
        rownames(m.beta.f) <- Estim.Opt$str.fixed.par
        
        ls.X.f <- lapply(ls.X, function(m.x){
            ##  IND*CT x NVAR
            m.x[, Estim.Opt$str.fixed.par, drop = FALSE]
        })
        
        ls.utility.f <- lapply(ls.X.f, function(m.x){
            m.u <- tcrossprod(m.x, t(v.beta.f))
            ##  IND*CT x CLASSES
            return(m.u)
        })
        
        ## IND*CT x DRAWS
        ls.utility <- mapply(function(m.x, m.y){m.x + m.y},
                             ls.utility, ls.utility.f, SIMPLIFY = FALSE)
        
        ########################################################################
        rm(ls.X.f, ls.utility.f)
        ########################################################################
    }
    
    ############################################################################
    ##  Check if we are estimating relative scale parameters
    ############################################################################
    if(Estim.Opt$b.relative.scale){
        v.lambda <- c(v.param[ls.str.par.names[["str.scale"]]], 1L)
        ##  IND*CT
        v.scale <- crossprod(t(m.R), v.lambda) / 1L
        v.scale <- as.vector(v.scale)
        
        ## IND*CT x DRAWS
        ls.utility <- lapply(ls.utility, function(m.x){
            return(m.x * v.scale)
        })
        
        ########################################################################
        rm(v.scale)
        ########################################################################
    }
    
    ############################################################################
    ##  Rescale utility
    ############################################################################
    if(Estim.Opt$b.rescale.utility){
        m.utility.max <- Reduce(pmax, ls.utility)
        m.utility.min <- Reduce(pmin, ls.utility)
        m.utility.mid <- (m.utility.max + m.utility.min) / 2L
        ls.utility <- lapply(ls.utility, function(m.x){m.x - m.utility.mid})
        
        ########################################################################
        rm(m.utility.max, m.utility.min, m.utility.mid)
        ########################################################################
    }
    
    ############################################################################
    ##  Calculating the probability of the chosen alternative
    ############################################################################
    ls.exp.utility <- lapply(ls.utility, function(m.x) exp(m.x))
    m.sum.utility <- Reduce("+", ls.exp.utility)
    ls.prob.alt <- lapply(ls.exp.utility, function(m.x) m.x /m.sum.utility)
    ls.prob.chosen <- mapply("*", ls.prob.alt, ls.Y, SIMPLIFY = FALSE)
    ##  IND*CT x DRAWS
    m.prob.chosen <- Reduce("+", ls.prob.chosen)
    m.prob.chosen <- matrix(m.prob.chosen, nrow = Estim.Opt$i.tasks)
    
    ############################################################################
    rm(ls.exp.utility, m.sum.utility, ls.prob.alt, ls.prob.chosen)
    ############################################################################
    
    ############################################################################
    ##  Check whether the data was complete
    ############################################################################
    if(!Estim.Opt$b.complete.data){
        if(any(is.nan(m.prob.chosen))){
            m.tmp <- as.logical(is.na(m.prob.chosen) - is.nan(m.prob.chosen))
            m.prob.chosen[m.tmp] <- 1L
            rm(m.tmp)
        } else{
            m.prob.chosen[is.na(m.prob.chosen)] <- 1L
        }
    }
    
    ############################################################################
    ##  Calculate the log likelihood value
    ############################################################################
    ##  Take the product over choice tasks - IND*DRAWS
    v.prob.sequence <- colProds(m.prob.chosen)
    
    ##  Average over draws - IND
    v.lik <- rowMeans2(matrix(v.prob.sequence, nrow = (nrow(ls.X[[1]])/Estim.Opt$i.tasks)))
    
    ##  Take the log
    v.log.lik <- log(v.lik)
    return(v.log.lik)
}
