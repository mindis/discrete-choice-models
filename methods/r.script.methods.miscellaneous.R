################################################################################
##  Name:         r.script.methods.miscellaneous.R
##  Created:      2018.08.02
##  Last edited:  2018.09.10
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for checking whether a function is byte-compiled
################################################################################
fnIsCompiled <- function(func){
    if(class(func) != "function") stop("only takes functions as arguments.")
    strLast2Lines <- taili(capture.output(func), 2)
    bByteCode <- any(grepl("bytecode:", strLast2Lines))
    return(bByteCode)
}

################################################################################
##  Function for searching best fitting model object in .GlobalEnv
################################################################################
fnFindBestFit <- function(EstimOpt){
    ##  Set some overall parameters
    iM <- EstimOpt$iM_run
    
    ##  Create a list containing all the model objects
    lsModel_max <- lapply(seq_len(iM), function(m){
        return(get(paste0("lsModel", m), envir = .GlobalEnv)$maximum)
    })
    
    ##  Put into a matrix and sort
    mOut <- matrix(NA, nrow = iM, ncol = 2L)
    m.tmp[, 1] <- seq_len(iM)
    m.tmp[, 2] <- Reduce(c, lsModel_max)
    
    ############################################################################
    rm(lsModel_max)
    ############################################################################
    
    ##  Sort and return
    mOut <- mOut[order(mOut[, 2L], decreasing = TRUE), ]
    return(mOut)
}