################################################################################
##  Name:         r.script.methods.data.R
##  Created:      2018.08.02
##  Last edited:  2018.09.10
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for generating complete data
################################################################################
fnCompleteData <- function(EstimOpt, mData){
    ##  Set some overall parameters
    iN <- EstimOpt$iN
    iJ <- EstimOpt$iJ
    iT <- EstimOpt$iT
    
    ##  Split the data into a list of matrices the length of EstimOpt$iN
    lsX <- vector(mode = "list", length = iN)
    lsX <- lapply(seq_len(iN), function(n){
        return(mData[mData[, EstimOpt$strID] == n, ])
    })
    
    ############################################################################
    ##  Function to repeat the data and create empty entries for missing CTs
    ############################################################################
    fn.complete.obs <- function(mX){
        ##  Check for possible coding errors
        if(nrow(mX) > (iJ * iT)){
            stop(paste("Number of rows per respondent", nrow(mX),
                       "is greater than",v(iJ * iT), ". Possible coding error.",
                       sep = " "))
        }
        
        ##  Check whether the data is complete for each respondent
        if(nrow(mX) == (iJ * iT)){
            return(mX)
        } else{
            ##  Create a NA matrix and fill in the data. NOTE: Doesn't handle CT order effects
            mX_tmp <- matrix(NA, nrow = (iJ * iT), ncol = ncol(mData))
            colnames(mX_tmp) <- colnames(mData)
            mX_tmp[, EstimOpt$strID] <- unique(mX[, EstimOpt$strID])
            mX_tmp[, EstimOpt$strCT] <- rep(seq_len(iT), each = iJ)
            mX_tmp[, EstimOpt$strALT] <- rep(seq_len(iJ), times = iT)
            mX_tmp[1L:nrow(mX), 4L:ncol(mX)] <- mX[, 4L:ncol(mX)]
            return(mX_tmp)
        }
    }
    
    ############################################################################
    ##  Run the function and return the complete data
    ############################################################################
    lsX <- lapply(lsX, function(mX){
        mOut <- fn.complete.obs(mX)
        return(mOut)
    })
    
    ############################################################################
    ##  Reduce, turn into matrix and return
    ############################################################################
    mData_complete <- as.matrix(Reduce(rbind, lsX))
    return(mData_complete)
}

################################################################################
##  Function for setting up the data
################################################################################
fnSetUpData <- function(EstimOpt){
    ##  Set some overall parameters
    iN <- EstimOpt$iN
    iJ <- EstimOpt$iJ
    iT <- EstimOpt$iT
    iR <- EstimOpt$iR
    
    ##  Read in the data -- only set up to work with RDS
    dfData <- readRDS(EstimOpt$strData)
    
    ##  Add a constant to the dataset if it does not exist
    if(!("const" %in% colnames(dfData))){
        const <- rep(1L, nrow(dfData))
        dfData <- cbind(dfData, const)
        rm(const)
    }
    
    ##  Subset the data to include only the variables entering the model
    strVars <- unique(c(EstimOpt$strID,
                        EstimOpt$strCT,
                        EstimOpt$strALT,
                        EstimOpt$strChoice,
                        EstimOpt$strP_fixed,
                        names(EstimOpt$lsP_rand),
                        unique(unlist(EstimOpt$lsP_het)),
                        EstimOpt$strP_class,
                        EstimOpt$strP_scale))
    
    dfData <- dfData[, strVars]
    
    ##  Turn into matrix for faster multiplication later
    mData <- as.matrix(dfData)
    colnames(mData) <- strVars
    
    ##  Check whether the data is complete and correct if not correct
    if(!EstimOpt$bCompleteData) mData <- fnCompleteData(EstimOpt, mData)
    
    ##  Split the data to be passed to the workers using a list of ids
    lsID <- split(mData[, EstimOpt$strID],
                  sort(mData[, EstimOpt$strID] %% EstimOpt$iCores))
    
    ##  IND*CT*ALT x NVAR -- adjusted for # of cores -- elements are matrices
    lsData <- lapply(lsID, function(vID){
        mX <- mData[mData[, EstimOpt$strID] %in% vID, ]
    })
    
    ############################################################################
    rm(dfData, strVars)
    invisible(gc(verbose = F))
    ############################################################################
    
    ############################################################################
    ##  Create an index for splitting the data to avoid a for-loop in the LL
    ############################################################################
    # lsData_index <- lapply(lsData, function(mX){
    #     lsX_tmp <- split(seq_len(nrow(mX[mX[, EstimOpt$strALT] == 1L, ])),
    #                     rep(1L:(nrow(mX)/(iT * iJ)), each = iT))
    # })
    # 
    # ##  Save the list to global environment
    # ##  IND*CT
    # assign("lsData_index", lsData_index, envir = .GlobalEnv)
    
    ############################################################################
    ##  Create a list of lists containing the attribute data - IND*CT x NVAR
    ############################################################################
    lsData_tmp <- lapply(lsData, function(mX){
        lsX_tmp <- vector(mode = "list", iJ)
        for(j in seq_len(iJ)){
            ##  IND*CT x NVAR
            lsX_tmp[[j]] <- mX[mX[, EstimOpt$strALT] == j,
                               c(EstimOpt$strP_fixed, names(EstimOpt$lsP_rand))]
        }
        
        ##  Name the alternatives
        names(lsX_tmp) <- paste0("alt", seq_len(iJ))
        
        ##  Check that the list elements are matrices
        lsX_tmp <- lapply(lsX_tmp, function(mX) {
            if(!is.matrix(mX)) mX <- as.matrix(mX)
            return(mX)
        })
        
        return(lsX_tmp)
    })
    
    ##  Save the list to global environment
    ##  IND*CT x NVAR 
    assign("lsDataX", lsData_tmp, envir = .GlobalEnv)
    
    ############################################################################
    ##  Create a list of lists containing the response data - IND*CT
    ############################################################################
    lsData_tmp <- lapply(lsData, function(mX){
        lsY_tmp <- vector(mode = "list", iJ)
        for(j in seq_len(iJ)){
            ##  IND*CT
            lsY_tmp[[j]] <- mX[mX[, EstimOpt$strALT] == j, EstimOpt$strChoice]
        }
        
        ##  Name the alternatives
        names(lsY_tmp) <- paste0("alt", seq_len(iJ))
        
        ##  Return the data IND*CT
        return(lsY_tmp)
    })
    
    ##  Save the list to global environment
    ##  IND*CT x CHOICE 
    assign("lsDataY", lsData_tmp, envir = .GlobalEnv)
    
    ############################################################################
    ##  Create a list of matrices containing interaction with means
    ##  IND*DRAWS x NVAR
    ############################################################################
    if(length(EstimOpt$lsP_het) > 0L){
        ##  Set some overall parameters
        strVars <- unique(unlist(EstimOpt$lsP_het))
        
        lsData_tmp <- lapply(lsData, function(mX){
            ##  IND*CT x NVAR
            mH_tmp <- mX[mX[, EstimOpt$strALT] == 1L, strVars, drop = F]
            
            ##  Create a list of row indices to use for splitting the data
            lsIndex <- split(seq_len(nrow(mH_tmp)), rep(1:(nrow(mH_tmp)/iT),
                                                        each = iT))
            
            ##  Average over CT and repeat equal to DRAWS
            lsH_tmp <- lapply(lsIndex, function(vX){
                mX_tmp <- mH_tmp[vX, , drop = F]
                ##  Added the na.rm = T to avoid an individual with missing CT data to have socio set to NA
                vX_tmp <- colMeans2(mX_tmp, na.rm = T)
                vX_tmp <- rep(vX_tmp, each = iR)
                mX <- matrix(vX_tmp, ncol = length(strVars))
                ##  DRAWS x NVAR
                return(mX)
            })
            
            ##  Reduce to matrix -- IND*DRAWS x NVAR
            mH_tmp <- Reduce(rbind, lsH_tmp)
            colnames(mH_tmp) <- strVars
            return(mH_tmp)
        })
        
        ##  Check that the list elements are matrices
        lsData_tmp <- lapply(lsData_tmp, function(mX) {
            if(!is.matrix(mX)) mX <- as.matrix(mX)
            return(mX)
        })
        
        ##  Save the list to global environment
        ##  IND*DRAWS x NVAR
        assign("lsDataH", lsData_tmp, envir = .GlobalEnv)
    }
    
    ############################################################################
    ##  Create a list of matrices containing relative scale variables
    ##  IND*CT x NVAR
    ############################################################################
    if(EstimOpt$bRelativeScale){
        ##  Set some overall parameters
        strVars <- EstimOpt$strP_scale
        
        lsData_tmp <- lapply(lsData, function(mX){
            ##  IND*CT x NVAR
            mX_tmp <- mX[mX[, EstimOpt$strALT] == 1L, strVars]
            colnames(mX_tmp) <- strVars
            return(mX_tmp)
        })
        
        ##  Check that the list elements are matrices
        lsData_tmp <- lapply(lsData_tmp, function(mX) {
            if(!is.matrix(mX)) mX <- as.matrix(mX)
            return(mX)
        })
        
        ##  Save the list to global environment
        ##  IND*CT x NVAR
        assign("lsDataR", lsData_tmp, envir = .GlobalEnv)
    }
    
    ############################################################################
    ##  Create a list of matrices for the latent class probabilities
    ##  IND*CT x NVAR
    ############################################################################
    if(EstimOpt$bLatentClass){
        ##  Set some overall parameters
        strVars <- EstimOpt$strP_class
        iQ <- EstimOpt$iQ
        
        ##  If we are estimating equality constraints we need to adjust iQ
        if(EstimOpt$bEqualityConstrained){
            if(EstimOpt$bDiscreteMixture){
                iQ <- length(EstimOpt$lsP_constrained)
            } else {
                if(!is.null(EstimOpt$mConstraints)){
                    iQ <- ncol(EstimOpt$mConstraints)
                } else {
                    iQ <- 2L ^ length(EstimOpt$lsP_constrained)
                }
            }
        }
        
        ##  Set up data for the class probability functions
        lsData_tmp <- lapply(lsData, function(mX){
            ##  Need a list the length of iQ
            lsC_tmp <- lapply(seq_len(iQ), function(q){
                ##  IND*CT x NVAR
                mC_tmp <- mX[mX[, EstimOpt$strALT] == 1L, strVars, drop = F]
                colnames(mC_tmp) <- strVars
                return(mC_tmp)
            })
            
            ##  Set the last vector (matrix) to zero, unless we are using a discrete mixture
            if(!EstimOpt$bDiscreteMixture){
                lsC_tmp[[iQ]] <- lsC_tmp[[iQ]] * 0L
            }
            
            ##  Check that the list elements are matrices
            lsC_tmp <- lapply(lsC_tmp, function(mX) {
                if(!is.matrix(mX)) mX <- as.matrix(mX)
                return(mX)
            })
            
            ##  IND*CT x NVAR
            return(lsC_tmp)
        })
        
        ##  Save the list to global environment
        ##  IND*CT x NVAR
        assign("lsDataC", lsData_tmp, envir = .GlobalEnv)
        
        ##  Now check if we are estimating a model with equality constraints
        if(EstimOpt$bEqualityConstrained){
            ##  Temporary vector of parameter names - NOTE: Doesn't handle combinations of fixed and random
            if(length(EstimOpt$lsP_rand) > 0L){
                strVars <- names(EstimOpt$lsP_rand)
            } else {
                strVars <- EstimOpt$strP_fixed
            }
            
            ##  Check whether we have a user supplied matrix
            if(is.null(EstimOpt$mConstraints)){
                ##  Expand to 2^k
                mC_tmp <- expand.grid(lapply(seq_along(EstimOpt$lsP_constrained),
                                             function(ik){return(c(0L, 1L))}))
            } else {
                mC_tmp <- t(EstimOpt$mConstraints)
                if(nrow(mC_tmp) < length(EstimOpt$lsP_constrained)){
                    stop("Check the dimensions of the supplied matrix of constraints! \n")
                }
            }
            
            ##  Assign the list of deltas to global
            lsC_tmp <- as.list(as.data.frame(t(mC_tmp)))
            assign("lsDelta", lsC_tmp, envir = .GlobalEnv)
            
            ##  If some attributes are dummy coded - repeat the relevant columns
            mC_tmp <- t(Reduce(cbind, lapply(seq_along(EstimOpt$lsP_constrained),
                                             function(k){
                                                 iK_tmp <- length(EstimOpt$lsP_constrained[[k]])
                                                 mC_tmp <- matrix(rep(mC_tmp[, k],
                                                                      times = iK_tmp),
                                                                  ncol = iK_tmp)
                                                 return(mC_tmp)
                                             })))
            
            rownames(mC_tmp) <- unlist(EstimOpt$lsP_constrained)
            
            ##  Sort based on rownames to match parameter vector
            # vSort <- rep(NA, nrow(mTemp))
            # for(k in seq_len(nrow(mTemp))){
            #     vSort[k] <- which(row.names(mTemp)[k] == strVars)
            # }
            # 
            # ##  NVAR x 2^k
            # mTemp <- mTemp[order(vSort), ]
            
            ##  Return as a list
            # lsTemp <- as.list(as.data.frame(mTemp))
            # lsTemp <- lapply(lsTemp, function(vX){
            #     names(vX) <- strVars
            #     return(vX)
            # })
            
            ##  NVAR x 2^k
            assign("mDeltaExpanded", mC_tmp, envir = .GlobalEnv)
        }
    }
    
    ############################################################################
    ##  Create a list containing the draws - IND*DRAWS x NVAR
    ############################################################################
    if(EstimOpt$bMakeDraws){
        ##  Create the matrix of draws and split to a list with elements equal
        ##  to the number of cores
        ##  IND*DRAWS x NVAR
        mDraws <- fnGenerateDraws(EstimOpt)
        lsCores <- vector(mode = "list", EstimOpt$iCores)
        lsCores <- lapply(seq_along(lsID), function(i){
            if(i == 1L) iE <<- 0L
            iS <- iE + (length(unique(lsID[[i]])) * iR)
            mDraws_tmp <- mDraws[(1L + iE):(iS), , drop = F]
            iE <<- iS
            return(mDraws_tmp)
        })
        
        ##  Check that the list elements are matrices
        lsCores <- lapply(lsCores, function(mX) {
            if(!is.matrix(mX)) mX <- as.matrix(mX)
            return(mX)
        })
        
        ##  Save the list to global environment
        ##  IND*DRAWS x NVAR
        assign("lsDraws", lsCores, envir = .GlobalEnv)
        
        ##  Create a list of draws indices
        # ls.tmp <- lapply(ls.tmp, function(mX){
        #     ls.tmp <- split(seq_len(nrow(mX)),
        #                     rep(1L:(nrow(mX)/iR),
        #                         each = iR))
        # })
        # 
        # ##  Save the list to global environment
        # ##  IND*DRAWS
        # assign("ls.draws.index", ls.tmp, envir = .GlobalEnv)
        
        ########################################################################
        rm(mDraws, lsCores)
        rm(iE, envir = .GlobalEnv)
        ########################################################################
    }
    
    ############################################################################
    rm(lsData, lsID, mData, lsData_tmp)
    invisible(gc(verbose = F))
    ############################################################################
    
}

