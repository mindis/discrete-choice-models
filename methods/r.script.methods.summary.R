################################################################################
##  Name:         r.script.methods.summary.R
##  Created:      2018.08.02
##  Last edited:  2018.09.10
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Create a function to calculate and print summary statistics
################################################################################
fnModelSummary <- function(lsM){
    ##  Set this option to avoid sceintific notation in printing
    # options(scipen = 999)
    
    ##  Set some overall parameters
    iN <- EstimOpt$iN
    iN_obs <- EstimOpt$iN_obs
    iJ <- EstimOpt$iJ
    iT <- EstimOpt$iT
    iM_search <- EstimOpt$iM_search
    iM_run <- EstimOpt$iM_run
    dLL_zero <- lsM$dLL_zero
    dLL <- lsM$maximum
    iK <- length(lsM$estimate)
    
    ############################################################################
    ##  Print general convergence information
    ############################################################################
    cat("###################################################################\n")
    cat("Model name: ", EstimOpt$strModelName, "\n")
    cat("Model message: ", lsM$message, "\n")
    cat("Convergence criteria: ", lsM$dConvergenceCriteria, "\n")
    cat("Estimating model ", lsM$m, "out of ",
        EstimOpt$iM_run, ". (", lsM$iCount_failed,
        "model(s) have failed).\n\n")
    
    cat("###################################################################\n")
    cat("Starting values (NOTE: Padded with zeros)\n")
    v_tmp <- c(lsM$vP_start, rep(0, ((5L * ceiling(iK / 5L)) - iK)))
    print(matrix(v_tmp, ncol = 5L, byrow = T))
    cat("\n")
    
    cat("###################################################################\n")
    cat("Final values (NOTE: Padded with zeros)\n")
    v_tmp <- c(lsM$estimate, rep(0, ((5L * ceiling(iK / 5L)) - iK)))
    print(matrix(v_tmp, ncol = 5L, byrow = T))
    cat("\n")
    
    cat("###################################################################\n")
    cat("Gradient at convergence (NOTE: Padded with zeros)\n")
    v_tmp <- c(lsM$gradient, rep(0, ((5L * ceiling(iK / 5L)) - iK)))
    print(matrix(v_tmp, ncol = 5L, byrow = T))
    cat("\n")
    
    ############################################################################
    ##  Create a function for printing some summary statistics
    ############################################################################
    fnSummaryStats <- function(dLL, dLL_zero, iK, iN_obs){
        vStats <- c(dLL,
                    dLL_zero,
                    iN_obs,
                    iK,
                    (1L - ((dLL - iK) / (dLL_zero))), 
                    ((-2L * dLL) + (2L * iK) ), 
                    ((-2L * dLL) + (3L * iK)), 
                    ((-2L * dLL) + (iK * (log(iN_obs) + 1L))), 
                    ((-2L * dLL) + (iK * (log((iN_obs + 2L) / 24L) + 1L ))), 
                    ((-2L * dLL) + (2L * iK)  +
                         (((2L * (iK + 1L)) * (iK + 2L))/(iN_obs - iK - 2L))),
                    ((-2L * dLL) + (iK * log(iN_obs))), 
                    ((-2L * dLL) + (iK * (log((iN_obs + 2L) / 24L)))), 
                    ((-2L * dLL) + (iK * (log(iN_obs) - log(2L * pi)))), 
                    ((-2L * dLL) + (2L * (iK * (log(log(iN_obs)))))))
        mStats_tmp <- matrix(vStats, nrow = 14L, ncol = 1L,
                             dimnames = list(c("LL", "LL(0)", "N", "K", "Adj. Rho^2",
                                               "AIC", "AIC3", "CAIC", "CAIC*",
                                               "HT-AIC/AICc", "BIC", "BIC*", "DBIC",
                                               "HQIC"), c("Value")))
        
        ##  Return the matrix
        return(mStats_tmp)
    }
    
    cat("###################################################################")
    cat("\n")
    cat("Model details: \n")
    mStats <- fnSummaryStats(dLL, dLL_zero, iK, iN_obs)
    print(mStats)
    cat("\n")
    
    ############################################################################
    ##  Calculate the variance covariance matrix and construct table of outputs
    ############################################################################
    ##  Check that we don't have any missing values in the Hessian matrix
    if(!any(is.nan(lsM$hessian))){
        mVCOV <- MASS::ginv(-lsM$hessian)
        
        ##  Check whether we are calculating robust standard errors
        if(EstimOpt$bRobustVCOV){
            mBread <- mVCOV * iN
            mBread[is.na(mBread)] <- 0L
            
            mMeat <- crossprod(lsM$gradientObs) / iN
            if(EstimOpt$bAdjustedRobustVCOV){
                mMeat <- mMeat * (iN / (iN - iK))
            }
            mMeat[is.na(mMeat)] <- 0L
            
            mVCOV <- (mBread %*% mMeat %*% mBread) / iN
        }
        
        ##  Matrix of outputs
        mOut <- matrix(NA, nrow = iK, ncol = 4L)
        mOut[, 1L] <- lsM$estimate
        mOut[, 2L] <- sqrt(diag(mVCOV))
        mOut[, 3L] <- mOut[, 1L] / mOut[, 2L]
        mOut[, 4L] <- 2L * pt(-abs(mOut[, 3L]), df = iN_obs)
        mOut <- round(mOut, 5)
        colnames(mOut) <- c("Est.", "Std. Err.", "T-stat", "P-value")
        rownames(mOut) <- names(lsM$estimate)
        
        cat("###################################################################")
        cat("\n")
        cat("Estimated parameters: \n")
        print(mOut)
        cat("\n")
        
        ########################################################################
        ##  Print extra information if we have a MIXL with correlations
        ########################################################################
        if(EstimOpt$bCorrelation == T){
            iK_f <- length(EstimOpt$strP_fixed)
            iK_r <- length(EstimOpt$lsP_rand)
            strVars <- names(EstimOpt$lsP_rand)
            iS <- (1L + iK_f + iK_r)
            iE <- (iK_f + ((2L * iK_r) + (iK_r * (iK_r - 1L) / 2L)))
            
            ##   Set up the lower Cholesky matrix
            mL <- matrix(0, iK_r, iK_r)
            mL[lower.tri(mL, diag = T)] <- lsM$estimate[iS:iE]
            
            ##   Create the correlation matrix
            mVCOV.param <- mL %*% t(mL)
            mD <- solve(diag(sqrt(diag(mVCOV.param))))
            mRho <- mD %*% mVCOV.param %*% mD
            
            ##   Combine the lower Cholesky and upper triangular correlation matrices
            m_tmp <- matrix(0, iK_r, iK_r)
            m_tmp[upper.tri(mRho, diag = F)] <- mRho[lower.tri(mRho, diag = F)]
            mCholCorr <- mL + m_tmp
            colnames(mCholCorr) <- strVars
            rownames(mCholCorr) <- strVars
            
            ##   Subset the vcov to match the lower Cholesky matrix
            mLambda <- mVCOV[iS:iE, iS:iE]
            
            ##   Calculate the SDs of the random parameters
            mSD <- fn.sd(mL, mLambda)
            colnames(mSD) <- c("Est.", "Std. Err.", "T-stat", "P-value")
            rownames(mSD) <- strVars
            
            ##   Calculate the standard errors of the correlations between the parameters
            mCorr <- fn.corr(mL, mLambda)
            colnames(mCorr) <-c("Est.", "Std. Err.", "T-stat", "P-value")
            strNames <- NULL
            
            for(i in 1:(length(EstimOpt$lsP_rand) - 1)){
                strNames <- c(strNames, paste(strVars[i],
                                              strVars[(i + 1):iK_r], sep = "."))
            }
            
            rownames(mCorr) <- strNames
            
            cat("###################################################################")
            cat("\n")
            cat("Standard deviations of the random parameters \n")
            print(mSD)
            cat("\n")
            
            cat("###################################################################")
            cat("\n")
            cat("Lower triangular Cholesky matrix and upper triangular correlation matrix \n")
            print(mCholCorr)
            cat("\n")
            
            cat("###################################################################")
            cat("\n")
            cat("Correlations between the random parameters \n")
            print(mCorr)
            cat("\n")
        }
        
        ########################################################################
        ##  Print the matrix of restrictions if we are estimating the ECLC models
        ########################################################################
        if(EstimOpt$bEqualityConstrained){
            strNames <- paste0("class.", seq_len(ncol(mDeltaExpanded)))
            m_tmp <- mDeltaExpanded
            colnames(m_tmp) <- strNames
            cat("###################################################################")
            cat("\n")
            cat("The class restrictions: \n")
            print(m_tmp)
            cat("\n")
        }
    } else {
        cat("###################################################################")
        cat("\n")
        cat("The Hessian contains NA. Moving to next model..")
        cat("\n")
    }
    
    
    tmp = fnTime(lsM$time)
    cat("###################################################################")
    cat("\n")
    cat("Model was estimated on", ifelse(EstimOpt$bParallel == T, EstimOpt$iCores, "1"), "core(s).\n")
    cat("Model completion time:", tmp[[1]], "days", tmp[[2]],"hours", tmp[[3]], "minutes and", tmp[[4]], "seconds \n")
    if(EstimOpt$bMakeDraws){
        cat("Estimated", length(EstimOpt$lsP_rand), "random variables using", EstimOpt$iR, ifelse(EstimOpt$bScramble == T, "scrambled", ""), EstimOpt$strDrawType, "draws.\n")
    }
    cat("The model is estimated in", ifelse(EstimOpt$bWTP_space == T, "WTP space.\n", "preference space.\n"))
    cat("The results are shown with", ifelse(EstimOpt$AadjustedRobustVCOV == T & EstimOpt$bRobustVCOV == T, "adjusted", ""), ifelse(EstimOpt$bRobustVCOV == T, "robust", "normal"), "standard errors.\n")
    cat("The model used the following seed: ", lsM$iSeed, "\n")
    cat(iN, "respondents made", EstimOpt$iN_obs, "choices.\n")
    cat("###################################################################")
    cat("\n")
}

################################################################################
##  Function for calculating time spent running the model
################################################################################
fnTime = function(t) {
    if(class(t) != "proc_time") stop("Input must be of class 'proc_time'")
    
    iTime <- as.numeric(t[3])
    
    iDays    <- floor( iTime / 86400)
    iHours   <- floor((iTime - (iDays * 86400)) / 3600)
    iMins <- floor((iTime - (iDays * 86400) - (iHours * 3600)) / 60)
    iSecs <- floor( iTime - (iDays * 86400) - (iHours * 3600) - (iMins * 60))
    
    return(list(iDays, iHours, iMins, iSecs))
}

################################################################################
##  Calculate the standard deviations and correlations of random parameters and
##  their corresponding standard erros
################################################################################

fn.1 = function(a, L){
    
    if(a > nrow(L)){stop(a," exceeds number of coefficients")}
    
    l = 1
    out = 0
    while(l < a + 1){
        out = out + L[a, l]^2
        l = l + 1
    }
    out
}

fn.2 = function(a, b, L){
    
    if( a >nrow(L)){stop(a," exceeds number of coefficients")}
    
    if(b > nrow(L)){stop(b," exceeds number of coefficients")}
    
    if(a > b){stop(a,">",b,": reverse order of coefficients")}
    
    l = 1
    out = 0
    
    while(l < a + 1){
        out = out + L[a, l] * L[b, l]
        l = l + 1
    }
    out
}

fn.3 = function(a, b, L){
    out = fn.2(a, b, L) / (sqrt(fn.1(a, L) * fn.1(b, L)))
    out
}

fn.1d = function(a, j, k, L){
    if(a > nrow(L)){stop(a, " exceeds number of coefficients")}
    if(k > j){stop(k, ">", j, ": k needs to be smaller than or equal to k")}
    
    out = 0
    if(a == j){
        out = 2 * L[a, k]
    }
    out
}

fn.1dv = function(a, L){
    
    K = sum(1:nrow(L))
    out = matrix(0, K, 1)
    
    j = 1
    n = 1
    while(j < nrow(L) + 1){
        k = 1
        while(k < j + 1){
            out[n, 1] = fn.1d(a, j, k, L)
            k = k + 1
            n = n + 1
        }
        j = j + 1
    }
    out
}

fn.2d = function(a, b, j, k, L){
    
    if(a > nrow(L)){stop(a, " exceeds number of coefficients")}
    if(b > nrow(L)){stop(b, " exceeds number of coefficients")}
    if(a > b){stop(a, ">", b, ": reverse order of coefficients")}
    if(a == b){stop(a, "=", b, ": use function with two separate coefficients coefficients")}
    if(k > j){stop(k, ">", j, ": k needs to be smaller than or equal to k")}
    
    l = 1
    out = 0
    
    if(a == j){
        out = L[b, k]
    }
    
    if(b == j){
        out = L[a, k]
    }
    out
}

fn.2dv = function(a, b, L){
    K = sum(1:nrow(L))
    out = matrix(0, K, 1)
    j = 1
    n = 1
    while(j < nrow(L) + 1){
        k = 1
        while(k < j + 1){
            out[n, 1] = fn.2d(a, b, j, k, L)
            k = k + 1
            n = n + 1
        }
        j = j + 1
    }
    out
}

fn.se.var = function(a, L, lambda){sqrt(t(fn.1dv(a, L)) %*% lambda %*% fn.1dv(a, L))}

fn.se.cor = function(a, b, L, lambda){sqrt(t(fn.2dv(a, b, L)) %*% lambda %*% fn.2dv(a, b, L))}

fn.3d = function(a, b, j, k, L){
    
    if(a > nrow(L)){stop(a, " exceeds number of coefficients")}
    if(b > nrow(L)){stop(b, " exceeds number of coefficients")}
    if(a > b){stop(a, ">", b, ": reverse order of coefficients")}
    if(a == b){stop(a, "=", b, ": use function with two separate coefficients coefficients")}
    if(k > j){stop(k, ">", j, ": k needs to be smaller than or equal to k")}
    
    l = 1
    out = (fn.2d(a, b, j, k, L) * sqrt(fn.1(a, L) * fn.1(b, L)) - (fn.2(a, b, L)) / (2 * (sqrt(fn.1(a, L) * fn.1(b, L)))) * (fn.1d(a, j, k, L) * fn.1(b, L) + fn.1d(b, j, k, L) * fn.1(a, L))) / (fn.1(a, L) * fn.1(b, L))
    out
}

fn.3dv = function(a, b, L){
    K = sum(1:nrow(L))
    out = matrix(0, K, 1)
    j = 1
    n = 1
    while(j < nrow(L) + 1){
        k = 1
        while(k < j + 1){
            out[n, 1] = fn.3d(a, b, j, k, L)
            k = k + 1
            n = n + 1
        }
        j = j + 1
    }
    out
}

fn.se.corr = function(a, b, L, lambda){sqrt(t(fn.3dv(a, b, L)) %*% lambda %*% fn.3dv(a, b, L))}

fn.4 = function(a, L){
    if(a > nrow(L)){stop(a, " exceeds number of coefficients")}
    l = 1
    out = 0
    while(l < a + 1){
        out = out + L[a, l]^2
        l = l + 1
    }
    out = sqrt(out)
    out
}

fn.4d = function(a, j, k, L){
    if(a > nrow(L)){stop(a, " exceeds number of coefficients")}
    if(k > j){stop(k, ">", j, ": k needs to be smaller than or equal to k")}
    
    out = 0
    if(a == j){
        out = 1 / (2 * fn.4(a, L)) * 2 * L[a, k]
    }
    out
}

fn.4dv = function(a, L){
    K = sum(1:nrow(L))
    out = matrix(0, K, 1)
    j = 1
    n = 1
    while(j < nrow(L) + 1){
        k = 1
        while(k < j + 1){
            out[n, 1] = fn.4d(a, j, k, L)
            k = k + 1
            n = n + 1
        }
        j = j + 1
    }
    out
}

fn.se.std = function(a, L, lambda){sqrt(t(fn.4dv(a, L)) %*% lambda %*% fn.4dv(a, L))}

fn.corr = function(L, lambda) {
    
    x = matrix(NA, sum(1:(nrow(L) - 1)), 4)
    i = 1
    n = 1
    
    while(i < nrow(L)){
        j = i + 1
        while(j < nrow(L) + 1) {
            
            y1 = fn.3(i, j, L)
            y2 = fn.se.corr(i, j, L, lambda)
            
            x[n, 1] = y1
            x[n, 2] = y2
            x[n, 3] = y1 / y2
            x[n, 4] = 2 * pt(-abs(y1 / y2), df = EstimOpt$iN_obs)
            
            j = j + 1
            n = n + 1
        }
        
        i = i + 1
    }
    
    return(x)
}

fn.sd = function(L, lambda){
    
    x = matrix(0, nrow(L), 4)
    i = 1
    
    while(i < nrow(L) + 1){
        
        y1 = fn.4(i, L)
        y2 = fn.se.std(i, L, lambda)
        
        x[i, 1] = y1
        x[i, 2] = y2
        x[i, 3] = y1 / y2
        x[i, 4] = 2 * pt(-abs(y1 / y2), df = EstimOpt$iN_obs)
        
        i = i + 1
        
    }
    return(x)
}




