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
    ##  Define some parameters
    i.Q <- length(ls.delta)
    
    ##  Check if we are calculating probabilities using a discrete mixture
    if(Estim.Opt$b.discrete.mixture){
        ##  This is the discrete mixture
        v.theta <- v.param[ls.str.par.names[["str.class"]]]
        
        ls.class.prob <- lapply(seq_along(ls.C), function(q){
            iK <- length(Estim.Opt$str.class.par)
            iS <- 1L + (iK * (q - 1L))
            iE <- iK * q
            ##  IND*CT
            return(as.vector(crossprod(t(ls.C[[q]]), v.theta[iS:iE])))
        })
        
        ##  Calculate the probability
        ls.class.prob <- lapply(ls.class.prob, function(vT){
            vPr <- 1L / (1L + exp(-(vT)))
            # vPr[is.na(vPr)] <- 0L
            return(vPr)
        })
        
        ##  PR(ATT) x IND*CT
        m.class.prob <- Reduce(rbind, ls.class.prob)
        
        ##  Use the probs to create the mixing distribution
        ls.class.prob <- lapply(ls.delta, function(vC){
            mT <- (m.class.prob * vC) + ((1L - m.class.prob) * (1L - vC))
            vT <- colProds(mT)
            return(vT)
        })
        
        ## IND*CT x CLASSES
        m.class.prob <- Reduce(cbind, ls.class.prob)
        
        ##  CT x IND*CLASSES
        m.class.prob <- matrix(as.vector(m.class.prob), nrow = Estim.Opt$i.tasks)
        
        ########################################################################
        rm(v.theta, ls.class.prob)
        ########################################################################
        
    } else {
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
        
        ########################################################################
        rm(v.theta, ls.class.prob, v.class.prob.sum)
        ########################################################################
    }
    
    ############################################################################
    ##  Make the betas IND*DRAWS x NVAR
    ############################################################################
    v.beta.r <- v.param[ls.str.par.names[["str.mean"]]]
    v.sigma.r <- v.param[ls.str.par.names[["str.std"]]]
    m.beta <- fn.make.beta(Estim.Opt, v.beta.r, v.sigma.r, m.draws)
    
    ############################################################################
    rm(v.beta.r, v.sigma.r)
    ############################################################################
    
    ############################################################################
    ##  Interactions with the means m.H -- IND*DRAWS x NVAR
    ############################################################################
    if(length(Estim.Opt$ls.het.par) > 0){
        v.phi <- v.param[ls.str.par.names[["str.het"]]]
        for(i in seq_along(Estim.Opt$ls.het.par)){
            str.tmp <- names(Estim.Opt$ls.het.par[i])
            v.phi.tmp <- v.phi[grep(str.tmp, names(v.phi))]
            ##  Multiply data with parameters
            m.tmp <- crossprod(t(m.H[, Estim.Opt$ls.het.par[[i]]]), v.phi.tmp)
            if(Estim.Opt$ls.rand.par[[str.tmp]] %in% c("-ln", "ln")){
                m.beta[, str.tmp] <- m.beta[, str.tmp] * exp(m.tmp)
            } else{
                m.beta[, str.tmp] <- m.beta[, str.tmp] + m.tmp
            }
        }
    }
    
    ############################################################################
    ##  Check if we are calculating utility in WTP space -- does not work with 
    ##  fixed parameters
    ############################################################################
    if(Estim.Opt$b.wtp.space){
        v.cost <- m.beta[, Estim.Opt$str.cost]
        m.beta[, Estim.Opt$str.cost] <- 1L
        ##  IND*DRAWS x NVAR
        m.beta <- m.beta * v.cost
    }
    
    ############################################################################
    ##  Create a list of betas (length of iQ) and calculate the random part of utility
    ############################################################################
    ls.beta <- lapply(seq_len(i.Q), function(q){
        mT <- m.beta
        ##  NVAR x IND*DRAWS
        return(t(mT))
    })
    
    ##  Keep an eye on this reverse argument to deal with the transpose
    ls.delta.expanded <- lapply(as.list(as.data.frame(m.delta.expanded)),
                                function(vX) {
                                    return(rev(vX))
                                })
    
    ##  Multiply with the restriction.
    ls.beta <- mapply("*", ls.beta, ls.delta.expanded, SIMPLIFY = F)
    
    ##  Constrained and of dim IND*DRAWS x NVAR
    ls.beta <- lapply(ls.beta, function(mB) return(t(mB)))
    
    ##  Subset the data to get the variables with random parameters
    ls.X.r <- lapply(ls.X, function(m.x){
        ##  IND*CT x NVAR
        m.x[, names(Estim.Opt$ls.rand.par), drop = FALSE]
    })
    
    i.T <- Estim.Opt$i.tasks
    i.D <- Estim.Opt$i.draws
    i.ind <- nrow(ls.X.r[[1L]]) / i.T
    
    ############################################################################
    ##  Use lapply() to calculate for each constrained beta matrix
    ############################################################################
    ls.lik <- lapply(ls.beta, function(mB){
        ########################################################################
        ##  Calculate the random part of utility
        ########################################################################
        ls.utility <- lapply(seq_len(Estim.Opt$i.alts), function(i.x){
            matrix(NA, nrow(ls.X.r[[1L]]), Estim.Opt$i.draws)
        })
        
        ##  IND*CT x DRAWS
        for(i.i in 1L:i.ind){
            v.rows <- (1L + ((i.i - 1L) * i.T)):(i.i * i.T)
            for(i.j in 1L:Estim.Opt$i.alts){
                ls.utility[[i.j]][v.rows, ] <- tcrossprod(ls.X.r[[i.j]][v.rows, , drop = FALSE],
                                                          mB[(1L + ((i.i - 1L) * i.D)):(i.i * i.D), ])
            }
        }
        
        ############################################################################
        ##  Calculate the fixed part of utility
        ############################################################################
        if(length(Estim.Opt$str.fixed.par) > 0){
            v.beta.f <- v.param[ls.str.par.names[["str.fixed"]]]
            ls.X.f <- lapply(ls.X, function(m.x){
                ##  IND*CT x NVAR
                m.x[, Estim.Opt$str.fixed.par, drop = FALSE]
            })
            
            ls.utility.f <- lapply(ls.X.f, function(m.x){
                v.u <- as.vector(crossprod(t(m.x), v.beta.f))
                ##  IND*CT
                return(v.u)
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
        return(v.lik)
    })
    
    ##  IND x CLASS
    m.lik <- Reduce(cbind, ls.lik)
  
    ##  IND*CLASS 
    v.class.prob <- colMeans2(m.class.prob, na.rm = T)
    ##  IND x CLASS
    m.class.prob <- matrix(v.class.prob, ncol = i.Q)
    
    ##  Calculate the likelihood
    v.lik <- rowSums2(m.class.prob * m.lik)
    
    ##  Take the log
    v.log.lik <- log(v.lik)
    return(v.log.lik)
}
