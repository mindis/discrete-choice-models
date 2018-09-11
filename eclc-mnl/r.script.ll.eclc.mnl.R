################################################################################
##  Name:         r.script.ll.eclc.mnl.R
##  Created:      2018.08.10
##  Last edited:  2018.09.11
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

fnLogLik <- function(vP){
    ##  Set some overall parameters
    iN <- nrow(lsX[[1]]) / EstimOpt$iT
    iT <- EstimOpt$iT
    iJ <- EstimOpt$iJ
    iR <- EstimOpt$iR
    iQ <- length(lsDelta)
    
    ############################################################################
    ##  Class probability function
    ############################################################################
    iK_c <- length(EstimOpt$strP_class)
    
    if(EstimOpt$bDiscreteMixture){
        vP_theta <- vP[lsParNames[["strP_class"]]]
        
        lsPr_att <- lapply(seq_along(lsC), function(q){
            iS <- 1L + (iK_c * (q - 1L))
            iE <- iK_c * q
            ##  IND*CT
            as.vector(crossprod(t(lsC[[q]]), vP_theta[iS:iE]))
        })
        
        ##  Calculate the probability
        lsPr_att <- lapply(lsPr_att, function(vC){
           1 / (1 + exp(-(vC)))
        })
        
        ##  Pr_att x IND*CT
        mPr_att <- Reduce(rbind, lsPr_att)
        
        ##  Calculate the class probabilities
        lsPr_class <- lapply(lsDelta, function(vD){
            mD <- (mPr_att * vD) + ((1 - mPr_att) * (1 - vD))
            vPr <- colProds(mD)
            return(vPr)
        })
        
        ##  IND*CT x iQ
        mPr_class <- Reduce(cbind, lsPr_class)
        mPr_class <- matrix(as.vector(mPr_class), nrow = iT)
        ##  IND x iQ - na.rm = T handles missing CTs
        mPr_class <- matrix(colMeans2(mPr_class, na.rm = T), ncol = iQ)
        
        ############################################################################
        rm(vP_theta, lsPr_att, lsPr_class, mPr_att)
        ############################################################################
        
    } else{
        ##  Zeros are added for the normalizing class
        vP_theta <- c(vP[lsParNames[["strP_class"]]], rep(0L, iK_c))
        
        lsPr_class <- lapply(seq_along(lsC), function(q){
            iS <- 1L + (iK_c * (q - 1L))
            iE <- iK_c * q
            ##  IND*CT
            exp(crossprod(t(lsC[[q]]), vP_theta[iS:iE]))
        })
        
        ##  Calculate the probability - Vector IND*CT
        vPr_class_sum <- Reduce("+", lsPr_class)
        lsPr_class <- lapply(lsPr_class, function(vC){
            vPr <- vC / vPr_class_sum
            # vPr[is.na(vPr)] <- 0
            as.vector(vPr)
        })
        
        ##  IND*CT x iQ
        mPr_class <- Reduce(cbind, lsPr_class)
        mPr_class <- matrix(as.vector(mPr_class), nrow = iT)
        ##  IND x iQ - na.rm = T handles missing CTs
        mPr_class <- matrix(colMeans2(mPr_class, na.rm = T), ncol = iQ)
        
        ############################################################################
        rm(vP_theta, lsPr_class, vPr_class_sum)
        ############################################################################
    }

    ############################################################################
    ##  Set up the matrix of betas
    ############################################################################
    mP_fixed <- matrix(rep(vP[lsParNames[["strP_fixed"]]], times = iQ),
                       ncol = iQ)
    rownames(mP_fixed) <- EstimOpt$strP_fixed
    
    ##  Impose the equality constraint
    mP_fixed <- mP_fixed * mDeltaExpanded
    
    ############################################################################
    ##  Check if we are calculating utility in WTP space 
    ############################################################################
    if(EstimOpt$bWTP_space){
        iC <- grep(EstimOpt$strP_cost, EstimOpt$strP_fixed)
        vP_cost <- mP_fixed[iC, ]
        mP_fixed[iC, ] <- 1L
        mP_fixed <- t(t(mP_fixed) * vP_cost)
    }    
    
    ############################################################################
    ##  Calculate the fixed part of utility
    ############################################################################
    lsU <- lapply(lsX, function(mX) tcrossprod(mX, t(mP_fixed)))
    
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
    lsPr_alt <- lapply(lsU_exp, function(mX) mX / mU_sum)
    lsPr_chosen <- mapply('*', lsPr_alt, lsY, SIMPLIFY = FALSE)
    mPr_chosen <- Reduce('+', lsPr_chosen)
    
    ##  CT x IND*iQ
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
    ##  Take the product over choice tasks - IND*iQ
    vPr_seq <- colProds(mPr_chosen)
    
    ##  Rearrange the LC matrix IND x CLASS
    mPr_seq <- matrix(vPr_seq, ncol = iQ)
    
    vLik <- rowSums2(mPr_class * mPr_seq)
    vLogLik <- log(vLik)
    return(vLogLik)
}