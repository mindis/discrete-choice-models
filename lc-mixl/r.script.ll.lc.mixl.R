################################################################################
##  Name:         r.script.ll.mixl.R
##  Created:      2018.08.02
##  Last edited:  2018.08.02
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

fn.log.lik <- function(v.param){
    ##  Get the parameters
    v.beta.r <- v.param[ls.str.par.names[["str.mean"]]]
    v.sigma.r <- v.param[ls.str.par.names[["str.std"]]]
    
    ############################################################################
    ##  Make the betas IND*DRAWS x NVAR
    ############################################################################
    m.beta <- fn.make.beta(Estim.Opt, v.beta.r, v.sigma.r, m.draws)
    
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
    
    ls.utility <- lapply(seq_len(Estim.Opt$i.alts), function(i.x){
        matrix(NA, nrow(ls.X.r[[1L]]), Estim.Opt$i.draws)
    })
    
    ##  IND*CT x DRAWS
    for(i.i in 1L:i.ind){
        v.rows <- (1L + ((i.i - 1L) * i.T)):(i.i * i.T)
        for(i.j in 1L:Estim.Opt$i.alts){
            ls.utility[[i.j]][v.rows, ] <- tcrossprod(ls.X.r[[i.j]][v.rows, , drop = FALSE],
                                                      m.beta[(1L + ((i.i - 1L) * i.D)):(i.i * i.D), ])
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
