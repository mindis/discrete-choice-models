################################################################################
##  Name:         r.script.methods.draws.R
##  Created:      2018.08.02
##  Last edited:  2018.09.10
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for generating MLHS draws
################################################################################
fnShuffle <- function(vInv){
    mOut <- vInv[rank(runif(length(vInv)))]
    return(mOut)
}

fnGenerateMLHS <- function(n, d, i){
    v_tmp <- seq(0, n - 1) / n
    mOut  <- matrix(0, n * i, d)
    
    j <- 1
    k <- 1
    
    while(j < i + 1) {
        k <- 1
        while(k < d + 1) {
            mOut[(1 + n * (j - 1)):(n * j), k] <- fnShuffle(v_tmp + runif(1)/n)
            k <- k + 1
        }
        j <- j + 1
    }
    return(mOut)
}

################################################################################
##  Function for generating scrambled Halton sequences
################################################################################
fnStart <- function(d, dMaxDigit){
    mCount <- matrix(0, d, dMaxDigit)
    mCount[, 1] <- rep(1, d)
    vDigit <- rep(1, d)
    
    return(list(mCount, vDigit))
}

##  Function for expanding the sequence of integers (Equation 1 in Bhat(2003))
fnDigitize <- function(d, vPrimes, mCount, vDigit){
    m <- 1L
    x <- NULL
    
    while(m <= d){
        l <- 1L
        r <- 1L
        while(r == 1L){
            x <- mCount[m, l] != (vPrimes[m] - 1L)
            r <- r - x
            mCount[m, l] <- (mCount[m, l] + 1L) * (x == 1L)
            vDigit[m] <- ((l - 1L) == vDigit[m]) + vDigit[m]
            l <- l + 1L
        }
        m <- m + 1L
    }
    return(list(mCount, vDigit))
}

##  Function for computing the radical inverse (Equation 2 in Bhat(2003))
fnRadInverse <- function(d, vPrimes, mCount, vDigit, mPermTable){
    m <- 1L
    mG <- matrix(0L, 1L, d)
    
    while(m <= d){
        l <- 1L
        dP <- vPrimes[m]
        while(l <= vDigit[m]){
            mG[m] <- (mPermTable[m, (mCount[m, l] + 1L)] / dP) + (mG[m])
            dP <- dP * vPrimes[m]
            l <- l + 1L
        }
        m <- m + 1L
    }
    return(mG)
}

##  Code for generating the Halton sequence
fnGenerateHalton <- function(n, d){
    if(d > 16L) stop("Cannot scramble Halton sequences beyond dimension 16.")
    
    dMaxDigit <- 50L
    vPrimes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53)
    mPermTable <- readRDS("../methods/m.xbrat.rds")
    
    mH <- matrix(0, n, d)
    vPrimes <- vPrimes[1:d]
    
    ##  Get values for the first number in the halton sequence
    ls.tmp <- fnStart(d, dMaxDigit)
    mCount <- ls.tmp[[1]]
    vDigit <- ls.tmp[[2]]
    
    mH[1L, ] <- fnRadInverse(d, vPrimes, mCount, vDigit, mPermTable)
    
    ##  For the rest of the numbers
    j <- 2L
    while(j <= n){
        ls.tmp <- fnDigitize(d, vPrimes, mCount, vDigit)
        mCount <- ls.tmp[[1]]
        vDigit <- ls.tmp[[2]]
        
        mH[j, ] <- fnRadInverse(d, vPrimes, mCount, vDigit, mPermTable)
        j <- j + 1L
    }
    return(mH)
}

################################################################################
##  Function for generating random draws
################################################################################
fnGenerateDraws <- function(EstimOpt){
    ##  Define some overall parameters
    iK <- length(EstimOpt$lsP_rand)
    iR <- EstimOpt$iR
    iN <- EstimOpt$iN
    iD <- EstimOpt$iD
    
    ##  If we are estimating LC-MIXL
    if(EstimOpt$bLatentClass && !EstimOpt$bEqualityConstrained){
        iK <- length(EstimOpt$lsP_rand) * EstimOpt$iQ
    }
    
    ##  Check that the specified type of draws are reckognized by the code
    if(!EstimOpt$strDrawType %in% c("pseudo", "halton", "sobol", "mlhs")){
        stop("Could not recognize type of draws. Please use one of the following: \n
             'pseudo', 'halton', 'sobol' or 'mlhs'.\n")
    }
    
    ##  Check that randtoolbox is loaded
    if(!("randtoolbox" %in% .packages())) require(randtoolbox)

    ##  Pseudo-random draws
    if(EstimOpt$strDrawType == "pseudo"){
        mDraws <- matrix(runif(iN * iR * iK), nrows = (iN * iR), ncol = iK)
    }
    
    ##  Halton draws
    if(EstimOpt$strDrawType == "halton"){
        if(EstimOpt$bScramble){
            mDraws <- fnGenerateHalton(((iN * iR) + iD), iK)
        } else {
            mDraws <- halton(((iN * iR) + iD), dim = iK)
        }
    }
    
    ##  Sobol draws
    if(EstimOpt$strDrawType == "sobol"){
        if(EstimOpt$bScramble){
            mDraws <- sobol((iN * iR), dim = iK, scrambling = EstimOpt$iScrambleType)
        } else {
            mDraws <- sobol((iN * iR), dim = iK)
        }
    }
    
    ##  MLHS draws
    if(EstimOpt$strDrawType == "mlhs"){
        mDraws <- fnGenerateMLHS(iN, iK, iR)
    }
    
    ##  Return the matrix of draws
    ## IND*DRAWS x NVAR
    mDraws <- as.matrix(mDraws)
    if((EstimOpt$strDrawType %in% c("halton") && iD > 0L)){
        return(as.matrix(qnorm(mDraws[-(1:iD), ])))
    } else {
        return(as.matrix(qnorm(mDraws)))
    }
}