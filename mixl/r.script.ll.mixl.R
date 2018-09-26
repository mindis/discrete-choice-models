################################################################################
##  Name:         r.script.ll.mixl.R
##  Created:      2018.08.02
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
    
    ############################################################################
    ##  Get the fixed parameters -- used to sort WTP space 
    ############################################################################
    if(length(EstimOpt$strP_fixed) > 0){
        vP_fixed <- vP[lsParNames[["strP_fixed"]]]
    }
    
    ############################################################################
    ##  Make the betas IND*DRAWS x NVAR
    ############################################################################
    vP_mean <- vP[lsParNames[["strP_mean"]]]
    vP_sigma <- vP[lsParNames[["strP_std"]]]
    mBeta <- fnMakeBeta(EstimOpt, vP_mean, vP_sigma, mDraws)
    
    ############################################################################
    ##  Interactions with the means mH -- IND*DRAWS x NVAR
    ############################################################################
    if(length(EstimOpt$lsP_het) > 0){
        vP_phi <- vP[lsParNames[["strP_het"]]]
        for(i in seq_along(EstimOpt$lsP_het)){
            strNames_tmp <- names(EstimOpt$lsP_het[i])
            vP_phi_tmp <- vP_phi[grep(strNames_tmp, names(vP_phi))]
            vZ_tmp <- crossprod(t(mH[, EstimOpt$lsP_het[[i]]]), vP_phi_tmp)
            if(EstimOpt$lsP_rand[[strNames_tmp]] %in% c("-ln", "ln")){
                mBeta[, strNames_tmp] <- mBeta[, strNames_tmp] * exp(vZ_tmp)
            } else{
                mBeta[, strNames_tmp] <- mBeta[, strNames_tmp] + vZ_tmp
            }
        }
    }
    
    ############################################################################
    ##  Check if we are calculating utility in WTP space -- does not work with 
    ##  fixed parameters
    ############################################################################
    if(EstimOpt$bWTP_space){
        vB_cost <- mBeta[, EstimOpt$strP_cost]
        mBeta[, EstimOpt$strP_cost] <- 1L
        ##  IND*DRAWS x NVAR
        mBeta <- mBeta * vB_cost
        rm(vB_cost)
    }
    
    ############################################################################
    ##  Calculate the random part of utility
    ############################################################################
    ##  Subset the data to get the variables with random parameters
    lsX_r <- lapply(lsX, function(mX){
        ##  IND*CT x NVAR
        mX[, names(EstimOpt$lsP_rand), drop = F]
    })
    
    ##  Empty matrix IND*CT x DRAWS
    lsU <- lapply(seq_len(iJ), function(i.x){
        matrix(NA, nrow(lsX_r[[1L]]), iR)
    })
    
    ##  IND*CT x DRAWS
    for(n in 1L:iN){
        vRows <- (1L + ((n - 1L) * iT)):(n * iT)
        for(j in 1L:iJ){
            lsU[[j]][vRows, ] <- tcrossprod(lsX_r[[j]][vRows, , drop = F],
                                            mBeta[(1L + ((n - 1L) * iR)):(n * iR), ])
        }
    }
    
    ############################################################################
    ##  Calculate the fixed part of utility
    ############################################################################
    if(length(EstimOpt$strP_fixed) > 0){
        lsX_f <- lapply(lsX, function(mX){
            ##  IND*CT x NVAR
            mX[, EstimOpt$strP_fixed, drop = F]
        })
        
        lsU_f <- lapply(lsX_f, function(mX){
            vU <- as.vector(crossprod(t(mX), vP_fixed))
            ##  IND*CT
            return(vU)
        })
        
        ## IND*CT x DRAWS
        lsU <- mapply(function(mX, mY){mX + mY}, lsU, lsU_f, SIMPLIFY = F)
        
        ########################################################################
        rm(lsX_f, lsU_f)
        ########################################################################
    }
    
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
    lsPr_alt <- lapply(lsU_exp, function(mX) mX / mU_sum )
    lsPr_chosen <- mapply('*', lsPr_alt, lsY, SIMPLIFY = FALSE)
    mPr_chosen <- Reduce('+', lsPr_chosen)
    
    ##  CT x IND*DRAWS
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
    ##  Take the product over choice tasks - IND*DRAWS
    vPr_seq <- colProds(mPr_chosen)
    
    ##  Average over draws - IND
    vLik <- rowMeans2(matrix(vPr_seq, nrow = iN))
    
    ##  Take the log
    vLogLik <- log(vLik)
    return(vLogLik)
}
