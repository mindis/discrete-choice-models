################################################################################
##  Name:         r.script.ll.mnl.R
##  Created:      2018.08.02
##  Last edited:  2018.08.02
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

fn.log.lik <- function(v.param){
    ##  Get the parameters
    v.beta.f <- v.param[ls.str.par.names[["str.fixed"]]]
    
    ############################################################################
    ##  Check if we are calculating utility in WTP space 
    ############################################################################
    if(Estim.Opt$b.wtp.space){
        i.cost.pos <- grep(Estim.Opt$str.cost, names(v.beta.f))
        d.cost <- v.beta.f[i.cost.pos]
        v.beta.f[i.cost.pos] <- 1L
        v.beta.f <- v.beta.f * d.cost
    }
    
    ############################################################################
    ##  Calculate the fixed part of utility
    ############################################################################
    ls.utility <- lapply(ls.X, function(m.x) crossprod(t(m.x), v.beta.f))
    
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
        v[is.na(v)] <- 0; as.vector(v)
    })
    ls.prob.chosen <- mapply('*', ls.prob.alt, ls.Y, SIMPLIFY = FALSE)
    m.prob.chosen <- as.matrix(Reduce('+', ls.prob.chosen))
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