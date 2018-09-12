################################################################################
##  Name:         r.script.ll.lc.mixl.R
##  Created:      2018.08.03
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
    iQ <- EstimOpt$iQ
    
    ############################################################################
    ##  Class probability function
    ############################################################################
    iK_c <- length(EstimOpt$strP_class)
    
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
    
    ############################################################################
    ##  Set up the fixed parameters -- if any
    ############################################################################
    if(length(EstimOpt$strP_fixed) > 0){
        mP_fixed <- matrix(vP[lsParNames[["strP_fixed"]]], ncol = iQ)
        rownames(mP_fixed) <- EstimOpt$strP_fixed
    }

    ############################################################################
    ##  Make a list of length equal to iQ containing the betas IND*DRAWS x NVAR
    ############################################################################
    mP_mean <- matrix(vP[lsParNames[["strP_mean"]]], ncol = iQ)
    mP_sigma <- matrix(vP[lsParNames[["strP_std"]]], ncol = iQ)
    rownames(mP_mean) <- rownames(mP_sigma) <- paste0("beta.",
                                                      names(EstimOpt$lsP_rand))
    
    iK_r <- nrow(mP_mean)
    
    ############################################################################
    ##  Wrap everyting inside a lapply() - loop
    ############################################################################
    lsLik <- lapply(seq_len(iQ), function(q){
        ########################################################################
        ##  Make the betas IND*DRAWS x NVAR
        ########################################################################
        iS <- 1 + (q - 1) * iK_r
        iE <- q * iK_r
        mBeta <- fnMakeBeta(EstimOpt, mP_mean[, q], mP_sigma[, q], mDraws[, iS:iE])
        
        ########################################################################
        ##  Interactions with the means mH -- IND*DRAWS x NVAR
        ########################################################################
        if(length(EstimOpt$lsP_het) > 0L){
            if(EstimOpt$bClassSpecific){
                iK_tmp <- length(lsParNames[["strP_het"]]) / iQ
                iS <- 1 + (q - 1) * iK_tmp
                iE <- q * iK_tmp
                vP_phi <- vP[lsParNames[["strP_het"]][iS:iE]]
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
            } else{
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
        }
        
        ########################################################################
        ##  Check if we are calculating utility in WTP space -- does not work with 
        ##  fixed parameters
        ########################################################################
        if(EstimOpt$bWTP_space){
            vB_cost <- mBeta[, EstimOpt$strP_cost]
            mBeta[, EstimOpt$strP_cost] <- 1L
            ##  IND*DRAWS x NVAR
            mBeta <- mBeta * vB_cost
            rm(vB_cost)
        }
        
        ########################################################################
        ##  Calculate the random part of utility
        ########################################################################
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
        
        ########################################################################
        ##  Calculate the fixed part of utility
        ########################################################################
        if(length(EstimOpt$strP_fixed) > 0){
            vP_fixed <- mP_fixed[, q]
            
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
            
            ####################################################################
            rm(lsX_f, lsU_f)
            ####################################################################
        }
        
        ########################################################################
        ##  Check if we are estimating relative scale parameters
        ########################################################################
        if(EstimOpt$bRelativeScale){
            vP_lambda <- c(vP[lsParNames[["strP_scale"]]], 1L)
            ##  IND*CT
            vS <- crossprod(t(mR), vP_lambda) / 1L
            vS <- as.vector(vS)
            
            ## IND*CT 
            lsU <- lapply(lsU, function(mX){
                return(mX * vS)
            })
            
            ####################################################################
            rm(vS)
            ####################################################################
        }
        
        ########################################################################
        ##  Rescale utility
        ########################################################################
        if(EstimOpt$bRescaleUtility){
            mU_max <- Reduce(pmax, lsU)
            mU_min <- Reduce(pmin, lsU)
            mU_mid <- (mU_max + mU_min) / 2L
            lsU <- lapply(lsU, function(mX){mX - mU_mid})
            
            ####################################################################
            rm(mU_max, mU_min, mU_mid)
            ####################################################################
        }
        
        ########################################################################
        ##  Calculating the probability of the chosen alternative
        ########################################################################
        lsU_exp <- lapply(lsU, function(mX) exp(mX))
        mU_sum <- Reduce("+", lsU_exp)
        lsPr_alt <- lapply(lsU_exp, function(mX) mX / mU_sum )
        lsPr_chosen <- mapply('*', lsPr_alt, lsY, SIMPLIFY = FALSE)
        mPr_chosen <- Reduce('+', lsPr_chosen)
        
        ##  CT x IND*DRAWS
        mPr_chosen <- matrix(mPr_chosen, nrow = iT)
        
        ########################################################################
        rm(lsU_exp, mU_sum, lsPr_alt, lsPr_chosen)
        ########################################################################
        
        ########################################################################
        ##  Check whether the data was complete
        ########################################################################
        if(!EstimOpt$bCompleteData){
            if(any(is.nan(mPr_chosen))){
                mPr_tmp <- as.logical(is.na(mPr_chosen) - is.nan(mPr_chosen))
                mPr_chosen[mPr_tmp] <- 1L
                rm(mPr_tmp)
            } else{
                mPr_chosen[is.na(mPr_chosen)] <- 1L
            }
        }
        
        ########################################################################
        ##  Calculate the log likelihood value
        ########################################################################
        ##  Take the product over choice tasks - IND*DRAWS
        vPr_seq <- colProds(mPr_chosen)
        
        ##  Average over draws - IND
        vLik <- rowMeans2(matrix(vPr_seq, nrow = iN))
        return(vLik)
        
    })

    ############################################################################
    ##  Calculate the log likelihood value
    ############################################################################
    ##  IND x iQ
    mLik <- Reduce(cbind, lsLik)
    
    ##  Calculate the likelihood
    vLik <- rowSums2(mPr_class * mLik)
    
    ##  Take the log
    vLogLik <- log(vLik)
    return(vLogLik)
}
