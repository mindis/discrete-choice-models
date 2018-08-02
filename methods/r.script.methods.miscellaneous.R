################################################################################
##  Name:         r.script.methods.miscellaneous.R
##  Created:      2018.08.02
##  Last edited:  2018.08.02
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for checking whether a function is byte-compiled
################################################################################
fn.is.compiled <- function(func){
    if(class(func) != "function") stop("only takes functions as arguments.")
    str.last.2.lines <- taili(capture.output(func), 2)
    b.byte.code <- any(grepl("bytecode:", str.last.2.lines))
    return(b.byte.code)
}

################################################################################
##  Function for searching best fitting model object in .GlobalEnv
################################################################################
fn.find.best.fit <- function(Estim.Opt){
    ##  Create a list containing all the model objects
    ls.model.max <- lapply(seq_len(Estim.Opt$i.nr.of.models), function(i.x){
        return(get(paste0("ls.model.", i.x), envir = .GlobalEnv)$maximum)
    })
    
    ##  Put into a matrix and sort
    m.tmp <- matrix(NA, nrow = Estim.Opt$i.nr.of.models, ncol = 2L)
    m.tmp[, 1] <- seq_len(Estim.Opt$i.nr.of.models)
    m.tmp[, 2] <- Reduce(c, ls.model.max)
    
    ############################################################################
    rm(ls.model.max)
    ############################################################################
    
    ##  Sort and return
    m.tmp <- m.tmp[order(m.tmp[, 2], decreasing = TRUE), ]
    return(m.tmp)
}