################################################################################
##  Name:         r.script.ll.lc.R
##  Created:      2018.08.10
##  Last edited:  2018.08.10
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

fn.log.lik <- function(v.param){
    ############################################################################
    ##  Class probability function
    ############################################################################
    ##  Define some parameters
    i.Q <- length(ls.constraints)
    
    ##  Check if we are calculating probabilities using a discrete mixture
    if(Estim.Opt$b.discrete.mixture){
        ##  This is the discrete mixture
        v.theta <- v.param[ls.str.par.names[["str.class"]]]
        
        ls.class.prob <- lapply(seq_along(ls.C), function(q){
            iK <- length(Estim.Opt$str.class.par)
            iS <- 1L + (iK * (q - 1L))
            iE <- iK * q
            ##  IND*CT
            exp(crossprod(t(ls.C[[q]]), v.theta[iS:iE]))
        })
        
        ##  Calculate the probability
        ls.class.prob <- lapply(ls.class.prob, function(vExp){
            vPr <- vExp / (vExp + 1L)
            vPr[is.na(vPr)] <- 0L
            return(as.vector(vPr))
        })
        
        ##  PR(ATT) x IND*CT
        m.class.prob <- t(Reduce(cbind, ls.class.prob))
        
        ##  Use the probs to create the mixing distribution
        ls.class.prob <- lapply(ls.constraints, function(vC){
            mT <- m.class.prob * vC + (1L - m.class.prob) * (1L - vC)
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
    ##  Set up the matrix of betas
    ############################################################################
    m.beta.f <- matrix(rep(v.param[ls.str.par.names[["str.fixed"]]], times = i.Q),
                       ncol = i.Q)
    rownames(m.beta.f) <- Estim.Opt$str.fixed.par
    
    ##  Impose the equality constraint
    m.const <- Reduce(cbind, ls.constraints)
    m.beta.f <- m.beta.f * m.const
    
    ############################################################################
    ##  Check if we are calculating utility in WTP space 
    ############################################################################
    if(Estim.Opt$b.wtp.space){
        i.cost.pos <- grep(Estim.Opt$str.cost, Estim.Opt$str.fixed.par)
        v.cost <- m.beta.f[i.cost.pos, ]
        m.beta.f[i.cost.pos, ] <- 1L
        m.beta.f <- t(t(m.beta.f) * v.cost)
    }    
    
    ############################################################################
    ##  Calculate the fixed part of utility
    ############################################################################
    ls.utility <- lapply(ls.X, function(m.x) tcrossprod(m.x, t(m.beta.f)))
    
    ############################################################################
    ##  Check if we are estimating relative scale parameters
    ############################################################################
    if(Estim.Opt$b.relative.scale){
        v.lambda <- c(v.param[ls.str.par.names[["str.scale"]]], 1L)
        ##  IND*CT
        v.scale <- crossprod(t(m.R), v.lambda) / 1L
        v.scale <- as.vector(v.scale)
        
        ## IND*CT 
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
    m.sum.utility <- Reduce('+', ls.exp.utility)
    ls.prob.alt <- lapply(ls.exp.utility, function(m.x) {
        v <- m.x / m.sum.utility
    })
    ls.prob.chosen <- mapply('*', ls.prob.alt, ls.Y, SIMPLIFY = FALSE)
    m.prob.chosen <- Reduce('+', ls.prob.chosen)
    ##  CT x IND*CLASS
    m.prob.chosen <- matrix(as.vector(m.prob.chosen), nrow = Estim.Opt$i.tasks)
    
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
    ##  IND*CLASS
    v.prob.sequence <- colProds(m.prob.chosen)
    ##  IND*CLASS
    v.class.prob <- colMeans2(m.class.prob)
    
    ##  Rearrange the LC matrix IND x CLASS
    m.prob.sequence <- matrix(v.prob.sequence, ncol = i.Q)
    m.class.prob <- matrix(v.class.prob, ncol = i.Q)
    
    v.lik <- rowSums2(m.class.prob * m.prob.sequence)
    v.log.lik <- log(v.lik)
    return(v.log.lik)
}