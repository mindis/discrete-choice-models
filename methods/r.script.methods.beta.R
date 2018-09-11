################################################################################
##  Name:         r.script.methods.beta.R
##  Created:      2018.08.02
##  Last edited:  2018.09.10
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for generating beta names
################################################################################
fnMakeBetaNames <- function(EstimOpt){
    ##  Define some overall parameters
    iQ <- EstimOpt$iQ
    
    ##  Create empty placeholders for the strings of parameter names
    strP_fixed <- strP_mean <- strP_std <- strP_het <- strP_scale <- strP_class <- NULL
    
    ##  Fixed parameters
    if(length(EstimOpt$strP_fixed) > 0){
        iK_f <- length(EstimOpt$strP_fixed)
        strP_fixed <- paste0("beta.", EstimOpt$strP_fixed)
        
        ##  If we are estimating LC-MNL and don't have equality constraints imposed
        if(EstimOpt$bLatentClass && !EstimOpt$bEqualityConstrained){
            iQ <- EstimOpt$iQ
            strP_fixed <- paste("beta", rep(EstimOpt$strP_fixed, times = iQ),
                                rep(seq_len(iQ), each = iK_f), sep = ".")
        }
    }
    
    ##  Random parameters
    if(length(EstimOpt$lsP_rand) > 0){
        strP_tmp <- names(EstimOpt$lsP_rand)
        strP_mean <- paste0("beta.", strP_tmp)
        strP_std <- paste0("sigma.", strP_tmp)
        if(EstimOpt$bCorrelation){
            strP_std <- NULL
            iK_r <- length(EstimOpt$lsP_rand)
            for(k in 1L:iK_r){
                strP_std <- c(strP_std,
                              paste("sigma", strP_tmp[k], strP_tmp[k:iK_r],
                                    sep = "."))
            }
        }
        
        ##  If we are estimating LC - MIXL and don't have equality constraints imposed
        if(EstimOpt$bLatentClass && !EstimOpt$bEqualityConstrained){
            strP_tmp <- names(EstimOpt$lsP_rand)
            strP_mean <- paste("beta", rep(strP_tmp, times = iQ),
                               rep(seq_len(iQ), each = length(strP_tmp)),
                               sep = ".")
            strP_std <- paste("sigma", rep(strP_tmp, times = iQ),
                              rep(seq_len(iQ), each = length(strP_tmp)),
                              sep = ".")
            if(EstimOpt$bCorrelation){
                strP_std <- NULL
                iK_r <- length(EstimOpt$lsP_rand)
                for(i in 1L:iK_r){
                    strP_std <- c(strP_std,
                                  paste("sigma", strP_tmp[i], strP_tmp[i:iK_r],
                                        sep = "."))
                }
                strP_std <- paste(rep(strP_std, times = iQ),
                                  rep(seq_len(iQ), each = length(strP_std)),
                                  sep = ".")
            }
        }
    }
    
    ##  Interactions with random parameters
    if(length(EstimOpt$lsP_het) > 0 ){
        strP_tmp <- names(EstimOpt$lsP_het)
        strP_het <- unlist(lapply(strP_tmp, function(x){
            outer(x, EstimOpt$lsP_het[[x]], FUN = paste, sep = ".")
        }))
        strP_het <- paste0("phi.", strP_het)
        ##  If we are estimating LC-MIXL and don't have equality constraints imposed
        if(EstimOpt$bLatentClass && EstimOpt$bClassSpecific && !EstimOpt$bEqualityConstrained){
            strP_het <- paste(rep(strP_het, times = iQ),
                              rep(seq_len(iQ), each = length(strP_het)),
                              sep = ".")
        }
    }
    
    ##  Relative scale parameters
    if(EstimOpt$bRelativeScale){
        iK_s <- length(EstimOpt$strP_scale)
        strP_scale <- paste0("lambda.", EstimOpt$strP_scale[-c(iK_s)])
    }
    
    ## Class probability function parameters
    if(EstimOpt$bLatentClass){
        iK_c <- length(EstimOpt$strP_class)
        iQ <- EstimOpt$iQ
        
        ## Change the number of classes if we are estimating equality constraints
        if(EstimOpt$bEqualityConstrained){
            if(EstimOpt$bDiscreteMixture){
                ##  Add the +1 to avoid another if - statement
                iQ <- length(EstimOpt$lsP_constrained) + 1L
            } else {
                if(!is.null(EstimOpt$mConstraints)){
                    iQ <- ncol(EstimOpt$mConstraints)
                } else {
                    iQ <- 2L^length(EstimOpt$lsP_constrained)
                }
            }
        }
        
        strP_class <- paste("theta",
                            rep(EstimOpt$strP_class, times = (iQ - 1L)),
                            rep(seq_len((iQ - 1L)), each = iK_c),
                            sep = ".")
    }
    
    ##  Return the list of parameter names
    return(list(strP_fixed = strP_fixed,
                strP_mean = strP_mean,
                strP_std = strP_std,
                strP_het = strP_het,
                strP_scale = strP_scale,
                strP_class = strP_class))
}

################################################################################
##  Function for generating betas
################################################################################
fnMakeBeta <- function(EstimOpt, vP_beta, vP_sigma, mDraws){
    ##  Check that the specified distributions are valid
    if(any(EstimOpt$lsP_rand %in% c("n", "ln", "-ln",
                                    "u", "-cu", "t",
                                    "ct", "sb") == FALSE)){
        stop("Cannot recognize the specified distributions.")
    }
    
    ##  IND*DRAWS x NVAR
    colnames(mDraws) <- names(EstimOpt$lsP_rand)
    
    ############################################################################
    ##  Function to transform the distributions
    ############################################################################
    fnTransformDistributions <- function(k){
        ##  Get the position of the draws for the corresponding beta
        iP <- which(names(vP_beta) == paste0("beta.", names(EstimOpt$lsP_rand)[k]))
        
        ##  Negative log-normal
        if(EstimOpt$lsP_rand[[k]] == "-ln"){
            if(EstimOpt$bCorrelation){
                vP_tmp <- exp(mBeta[, iP]) * -1
            } else {
                vP_tmp <- exp(vP_beta[iP] + vP_sigma[iP] * mDraws[, iP, drop = FALSE]) * -1
            }
        }
        
        ##  Log-normal
        if(EstimOpt$lsP_rand[[k]] == "ln"){
            if(EstimOpt$bCorrelation){
                vP_tmp <- exp(mBeta[, iP])
            } else {
                vP_tmp <- exp(vP_beta[iP] + vP_sigma[iP] * mDraws[, iP, drop = FALSE])
            }
        }
        
        ##  Normal
        if(EstimOpt$lsP_rand[[k]] == "n"){
            if(EstimOpt$bCorrelation){
                vP_tmp <- mBeta[, iP]
            } else {
                vP_tmp <- vP_beta[iP] + vP_sigma[iP] * mDraws[, iP, drop = FALSE]
            }
        }
        
        ##  Uniform -1, 1
        if(EstimOpt$lsP_rand[[k]] == "u"){
            vDraws_tmp <- 2 * pnorm(mDraws[, iP, drop = FALSE]) - 1
            vP_tmp <- vP_beta[iP] + vP_sigma[iP] * vDraws_tmp
        }
        
        ##  Constrained negative uniform
        if(EstimOpt$lsP_rand[[k]] == "-cu"){
            vDraws_tmp <- 2 * pnorm(mDraws[, iP, drop = FALSE]) - 1
            vP_tmp <- exp(vP_beta[iP] + vP_sigma[iP] * vDraws_tmp) * -1
        }
        
        ##  Triangular
        if(EstimOpt$lsP_rand[[k]] == "t"){
            vDraws_tmp <- pnorm(mDraws[, iP, drop = FALSE])
            bDraws_tmp <- vDraws_tmp < 0.5
            vDraws_tmp <- bDraws_tmp * (sqrt(2 * vDraws_tmp) - 1) +
                !bDraws_tmp * (1 - sqrt(2 * (1 - vDraws_tmp)))
            vP_tmp <- vP_beta[iP] + vP_sigma[iP] * vDraws_tmp
        }
        
        ##  Constrained triangular
        if(EstimOpt$lsP_rand[[k]] == "ct"){
            vDraws_tmp <- pnorm(mDraws[, iP, drop = FALSE])
            bDraws_tmp <- vDraws_tmp < 0.5
            vDraws_tmp <- bDraws_tmp * (sqrt(2 * vDraws_tmp) - 1) +
                !bDraws_tmp * (1 - sqrt(2 * (1 - vDraws_tmp)))
            d.spread <- vP_beta[iP] * (1 / (1 + exp(vP_sigma[iP])))
            vP_tmp <- vP_beta[iP] + d.spread * vDraws_tmp
        }
        
        ##  Johnson SB -- Not set up to work properly
        # if(EstimOpt$lsP_rand[[k]] == "sb"){
        #     vP_tmp <- vP_beta[iP] + vP_sigma[iP] * mDraws[, iP, drop = FALSE]
        #     vP_tmp <- exp(vP_tmp) / (1 + exp(vP_tmp))
        # }
        
        ##  Return the transformed parameter
        return(vP_tmp)
    }
    
    ############################################################################
    ##  Check for correlation
    ############################################################################
    if(EstimOpt$bCorrelation){
        if(any(EstimOpt$lsP_rand %in% c("-ln", "ln", "n") == FALSE)){
            stop("Only distributions in the family of normals can be correlated")
        }
        
        ##  Set up the lower Cholesky matrix
        iK <- length(EstimOpt$lsP_rand)
        mL <- matrix(0L, iK, iK)
        mL[lower.tri(mL, diag = TRUE)] <- vP_sigma
        
        ##  Calculate the betas
        mSigma <- tcrossprod(mL, mDraws)
        
        ##  IND*DRAWS x NVAR
        mBeta <- vP_beta + mSigma
        mBeta <- t(mBeta)
    } else {
        ##  Transform the distributions
        lsBeta <- lapply(seq_along(EstimOpt$lsP_rand), function(k){
            vB <- fnTransformDistributions(k)
        })
        
        mBeta <- Reduce(cbind, lsBeta)
    }
    
    ## IND*DRAWS x NVAR
    if(!is.matrix(mBeta)) mBeta <- as.matrix(mBeta)
    colnames(mBeta) <- names(EstimOpt$lsP_rand)
    return(mBeta)
}
