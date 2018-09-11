################################################################################
##  Name:         r.script.ll.mnl.R
##  Created:      2018.08.02
##  Last edited:  2018.09.10
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

fnLogLik <- function(vP){
    ##  Set some overall parameters
    iN <- nrow(lsX[[1]]) / EstimOpt$iT
    iT <- EstimOpt$iT
    iJ <- EstimOpt$iJ
    
    ##  Get the parameters
    vP_fixed <- vP[lsParNames[["strP_fixed"]]]
    
    ############################################################################
    ##  Check if we are calculating utility in WTP space 
    ############################################################################
    if(EstimOpt$bWTP_space){
        iC <- grep(EstimOpt$strP_cost, names(vP_fixed))
        dP_tmp <- vP_fixed[iC]
        vP_fixed[iC] <- 1L
        vP_fixed <- vP_fixed * dP_tmp
    }
    
    ############################################################################
    ##  Calculate the fixed part of utility
    ############################################################################
    lsU <- lapply(lsX, function(mX) crossprod(t(mX), vP_fixed))
    
    ############################################################################
    ##  Check if we are estimating relative scale parameters
    ############################################################################
    if(EstimOpt$bRelativeScale){
        vP_lambda <- c(vP[lsParNames[["strP_scale"]]], 1L)
        ##  IND*CT
        vS <- crossprod(t(mR), vP_lambda) / 1L
        vS <- as.vector(vS)
        
        ## IND*CT 
        lsU <- lapply(lsU, function(mX){
            return(mX * vS)
        })
        
        ########################################################################
        rm(vS)
        ########################################################################
    }
    
    ############################################################################
    ##  Rescale utility
    ############################################################################
    if(EstimOpt$bRescaleUtility){
        mU_max <- Reduce(pmax, lsU)
        mU_min <- Reduce(pmin, lsU)
        mU_mid <- (mU_max + mU_min) / 2L
        lsU <- lapply(lsU, function(mX){mX - mU_mid})
        
        ########################################################################
        rm(mU_max, mU_min, mU_mid)
        ########################################################################
    }
    
    ############################################################################
    ##  Calculating the probability of the chosen alternative
    ############################################################################
    lsU_exp <- lapply(lsU, function(mX) exp(mX))
    mU_sum <- Reduce("+", lsU_exp)
    lsPr_alt <- lapply(lsU_exp, function(mX) {
        v <- mX / mU_sum
        v[is.na(v)] <- 0
        return(as.vector(v))
    })
    
    lsPr_chosen <- mapply('*', lsPr_alt, lsY, SIMPLIFY = FALSE)
    mPr_chosen <- Reduce('+', lsPr_chosen)
    
    ##  CT x IND
    mPr_chosen <- matrix(mPr_chosen, nrow = iT)
    
    ############################################################################
    rm(lsU_exp, mU_sum, lsPr_alt, lsPr_chosen)
    ############################################################################
    
    ############################################################################
    ##  Check whether the data was complete
    ############################################################################
    if(!EstimOpt$bCompleteData){
        if(any(is.nan(mPr_chosen))){
            mPr_tmp <- as.logical(is.na(mPr_chosen) - is.nan(mPr_chosen))
            mPr_chosen[mPr_tmp] <- 1L
            rm(mPr_tmp)
        } else{
            mPr_chosen[is.na(mPr_chosen)] <- 1L
        }
    }
    
    ############################################################################
    ##  Calculate the log likelihood value
    ############################################################################
    ##  Take the product over choice tasks - IND
    vPr_seq <- colProds(mPr_chosen)
    
    ##  Average over draws - IND -- for the MNL this does not matter
    vLik <- rowMeans2(matrix(vPr_seq, nrow = iN))
    
    ##  Take the log
    vLogLik <- log(vLik)
    return(vLogLik)
}