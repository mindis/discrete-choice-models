################################################################################
##  Name:         r.script.methods.data.R
##  Created:      2018.08.02
##  Last edited:  2018.08.02
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for generating complete data
################################################################################
fn.complete.data <- function(Estim.Opt, m.data){
    ##  Split the data into a list of matrices the length of Estim.Opt$i.ind
    ls.X <- vector(mode = "list", length = Estim.Opt$i.ind)
    ls.X <- lapply(seq_len(Estim.Opt$i.ind), function(i.x){
        return(m.data[m.data[, Estim.Opt$str.id] == i.x, ])
    })
    
    ##  Cretate indices and a string of variable names
    i.rows <- Estim.Opt$i.alts * Estim.Opt$i.tasks
    i.cols <-  ncol(m.data)
    str.names <- colnames(m.data)
    
    ############################################################################
    ##  Function to repeat the data and create empty entries for missing CTs
    ############################################################################
    fn.complete.obs <- function(m.x){
        ##  Check for possible coding errors
        if(nrow(m.x) > i.rows){
            stop(paste("Number of rows per respondent", x, "is greater than",
                       i.rows, ". Possible coding error.", sep = " "))
        }
        
        ##  Check whether the data is complete for each respondent
        if(nrow(m.x) == i.rows){
            return(m.x)
        } else{
            ##  Create an empty matrix of NA of correct dimensions and fill in
            ##  the data that we do have -- NA is handled in the LL function
            m.tmp <- matrix(NA, nrow = i.rows, ncol = i.cols)
            colnames(m.tmp) <- str.names
            m.tmp[, Estim.Opt$str.id] <- unique(m.x[, Estim.Opt$str.id])
            m.tmp[, Estim.Opt$str.ct] <- rep(seq_len(Estim.Opt$i.tasks),
                                             each = Estim.Opt$i.alts)
            m.tmp[, Estim.Opt$str.alt] <- rep(seq_len(Estim.Opt$i.alts),
                                              times = Estim.Opt$i.tasks)
            ##  The first three columns are always id, ct, alt
            ##  Note! Ordering effects of ct's are not considered
            m.tmp[1L:nrow(m.x), 4L:ncol(m.x)] <- m.x[, 4L:ncol(m.x)]
            return(m.tmp)
        }
    }
    
    ############################################################################
    ##  Run the function and return the complete data
    ############################################################################
    ls.X <- lapply(ls.X, function(m.x){
        m.out <- fn.complete.obs(m.x)
        return(m.out)
    })
    
    ############################################################################
    ##  Reduce, turn into matrix and return
    ############################################################################
    m.data.out <- as.matrix(Reduce(rbind, ls.X))
    return(m.data.out)
}

################################################################################
##  Function for setting up the data
################################################################################
fn.set.up.data <- function(Estim.Opt){
    ##  Read in the data -- Only set up to work with RDS
    df.data <- readRDS(Estim.Opt$str.data)
    
    ##  Add a constant to the dataset if it does not exist
    if(!("const" %in% colnames(df.data))){
        const <- rep(1L, nrow(df.data))
        df.data <- cbind(df.data, const)
        rm(const)
    }
    
    ##  Subset the data to include only the variables entering the model
    str.var.names <- unique(c(Estim.Opt$str.id,
                              Estim.Opt$str.ct,
                              Estim.Opt$str.alt,
                              Estim.Opt$str.choice,
                              Estim.Opt$str.fixed.par,
                              names(Estim.Opt$ls.rand.par),
                              unique(unlist(Estim.Opt$ls.het.par)),
                              Estim.Opt$str.class.par,
                              Estim.Opt$str.scale.par))
    
    df.data <- df.data[, str.var.names]
    
    ##  Turn into matrix for faster multiplication later
    m.data <- as.matrix(df.data)
    colnames(m.data) <- str.var.names
    
    ##  Check whether the data is complete and correct if not correct
    if(!Estim.Opt$b.complete.data) m.data <- fn.complete.data(Estim.Opt, m.data)
    
    ##  Split the data to be passed to the workers using a list of ids
    ls.id <- split(m.data[, Estim.Opt$str.id],
                   sort(m.data[, Estim.Opt$str.id] %% Estim.Opt$i.cores))
    
    ##  IND*CT*ALT x NVAR -- adjusted for # of cores -- elements are matrices
    ls.data <- lapply(ls.id, function(v.id){
        m.X <- m.data[m.data[, Estim.Opt$str.id] %in% v.id, ]
    })
    
    ############################################################################
    rm(df.data, str.var.names)
    invisible(gc(verbose = FALSE))
    ############################################################################
    
    ############################################################################
    ##  Create an index for splitting the data to avoid a for-loop in the LL
    ############################################################################
    ls.data.tmp <- lapply(ls.data, function(m.x){
        ls.tmp <- split(seq_len(nrow(m.x[m.x[, Estim.Opt$str.alt] == 1L, ])),
                        rep(1L:(nrow(m.x)/(Estim.Opt$i.tasks * Estim.Opt$i.alts)),
                                each = Estim.Opt$i.tasks))
    })
    
    ##  Save the list to global environment
    ##  IND*CT
    assign("ls.data.index", ls.data.tmp, envir = .GlobalEnv)
    
    ############################################################################
    ##  Create a list of lists containing the attribute data - IND*CT x NVAR
    ############################################################################
    ls.data.tmp <- lapply(ls.data, function(m.x){
        ls.tmp <- vector(mode = "list", Estim.Opt$i.alts)
        for(i.j in seq_len(Estim.Opt$i.alts)){
            ##  IND*CT x NVAR
            ls.tmp[[i.j]] <- m.x[m.x[, Estim.Opt$str.alt] == i.j,
                                 c(Estim.Opt$str.fixed.par,
                                   names(Estim.Opt$ls.rand.par))]
        }
        
        ##  Name the alternatives
        names(ls.tmp) <- paste0("alt", seq_len(Estim.Opt$i.alts))
        
        ##  Check that the list elements are matrices
        ls.tmp <- lapply(ls.tmp, function(m.x) {
            if(!is.matrix(m.x)) m.x <- as.matrix(m.x)
            return(m.x)
        })
        
        return(ls.tmp)
    })
    
    ##  Save the list to global environment
    ##  IND*CT x NVAR 
    assign("ls.data.X", ls.data.tmp, envir = .GlobalEnv)
    
    ############################################################################
    ##  Create a list of lists containing the response data - IND*CT
    ############################################################################
    ls.data.tmp <- lapply(ls.data, function(m.x){
        ls.tmp <- vector(mode = "list", Estim.Opt$i.alts)
        for(i.j in seq_len(Estim.Opt$i.alts)){
            ##  IND*CT
            ls.tmp[[i.j]] <- m.x[m.x[, Estim.Opt$str.alt] == i.j,
                                 Estim.Opt$str.choice]
        }
        
        ##  Name the alternatives
        names(ls.tmp) <- paste0("alt", seq_len(Estim.Opt$i.alts))
        
        ##  Check that the list elements are matrices
        # ls.tmp <- lapply(ls.tmp, function(m.x) {
        #     if(!is.matrix(m.x)) m.x <- as.matrix(m.x)
        #     return(m.x)
        # })
        
        ##  Return the data IND*CT
        return(ls.tmp)
    })
    
    ##  Save the list to global environment
    ##  IND*CT x CHOICE 
    assign("ls.data.Y", ls.data.tmp, envir = .GlobalEnv)
    
    ############################################################################
    ##  Create a list of matrices containing interaction with means
    ##  IND*DRAWS x NVAR
    ############################################################################
    if(length(Estim.Opt$ls.het.par) > 0L){
        str.col.names <- unique(unlist(Estim.Opt$ls.het.par))
        ls.data.tmp <- lapply(ls.data, function(m.x){
            ##  IND*CT x NVAR
            m.tmp <- m.x[m.x[, Estim.Opt$str.alt] == 1L, str.col.names,
                         drop = FALSE]
            
            ##  Create a list of row indices to use for splitting the data
            ls.index <- split(seq_len(nrow(m.tmp)),
                              rep(1:(nrow(m.tmp)/Estim.Opt$i.tasks),
                                  each = Estim.Opt$i.tasks))
            
            ##  Average over CT and repeat equal to DRAWS
            ls.tmp <- lapply(ls.index, function(v.x){
                m.x <- m.tmp[v.x, , drop = FALSE]
                v.tmp <- colMeans2(m.x)
                v.tmp <- rep(v.tmp, each = Estim.Opt$i.draws)
                m.tmp <- matrix(v.tmp, ncol = length(str.col.names))
                ##  DRAWS x NVAR
                return(m.tmp)
            })
            
            
            ##  Average over CT and repeat equal to DRAWS
            # ls.tmp <- lapply(ls.tmp, function(m.x){
            #     m.x <- matrix(m.x, ncol = length(str.col.names))
            #     if(!is.matrix(m.x)) m.x <- as.matrix(m.x)
            #     v.tmp <- colMeans2(m.x)
            #     m.tmp <- matrix(v.tmp, nrow = Estim.Opt$i.draws, byrow = TRUE)
            #     ##  DRAWS x NVAR
            #     return(m.tmp)
            # })
            
            ##  Reduce to matrix -- IND*DRAWS x NVAR
            m.tmp <- Reduce(rbind, ls.tmp)
            colnames(m.tmp) <- str.col.names
            
            return(m.tmp)
        })
        
        ##  Check that the list elements are matrices
        ls.data.tmp <- lapply(ls.data.tmp, function(m.x) {
            if(!is.matrix(m.x)) m.x <- as.matrix(m.x)
            return(m.x)
        })
        
        ##  Save the list to global environment
        ##  IND*DRAWS x NVAR
        assign("ls.data.H", ls.data.tmp, envir = .GlobalEnv)
    }
    
    ############################################################################
    ##  Create a list of matrices containing relative scale variables
    ##  IND*CT x NVAR
    ############################################################################
    if(Estim.Opt$b.relative.scale){
        str.col.names <- Estim.Opt$str.scale.par
        ls.data.tmp <- lapply(ls.data, function(m.x){
            ##  IND*CT x NVAR
            m.tmp <- m.x[m.x[, Estim.Opt$str.alt] == 1L, str.col.names]
            colnames(m.tmp) <- str.col.names
            return(m.tmp)
        })
        
        ##  Check that the list elements are matrices
        ls.data.tmp <- lapply(ls.data.tmp, function(m.x) {
            if(!is.matrix(m.x)) m.x <- as.matrix(m.x)
            return(m.x)
        })
        
        ##  Save the list to global environment
        ##  IND*CT x NVAR
        assign("ls.data.R", ls.data.tmp, envir = .GlobalEnv)
    }
    
    ############################################################################
    ##  Create a list of matrices for the latent class probabilities
    ##  IND*CT x NVAR
    ############################################################################
    if(Estim.Opt$b.latent.class){
        str.col.names <- Estim.Opt$str.class.par
        i.Q <- Estim.Opt$i.classes
        
        ls.data.tmp <- lapply(ls.data, function(m.x){
            ##  Need a list the length of i.Q
            ls.tmp <- lapply(seq_len(i.Q), function(i.q){
                ##  IND*CT x NVAR
                m.tmp <- m.x[m.x[, Estim.Opt$str.alt] == 1L, str.col.names,
                             drop = FALSE]
                colnames(m.tmp) <- str.col.names
                return(m.tmp)
            })
            
            ##  Set the last vector (matrix) to zero
            ls.tmp[[i.Q]] <- ls.tmp[[i.Q]] * 0L
            
            ##  Check that the list elements are matrices
            ls.tmp <- lapply(ls.tmp, function(m.x) {
                if(!is.matrix(m.x)) m.x <- as.matrix(m.x)
                return(m.x)
            })
            
            ##  IND*CT x NVAR
            return(ls.tmp)
        })
        
        ##  Save the list to global environment
        ##  IND*CT x NVAR
        assign("ls.data.C", ls.data.tmp, envir = .GlobalEnv)
        
        ##  Now check if we are estimating a model with equality constraints
        if(Estim.Opt$b.equality.constrained){
            ##  Temporary vector of parameter names
            strT <- ifelse(length(Estim.Opt$ls.rand.par) > 0L,
                           names(Estim.Opt$ls.rand.par),
                           Estim.Opt$str.fixed.par)
            
            ##  Check whether we have a user supplied matrix
            if(is.null(Estim.Opt$m.constraints)){
                ##  Expand to 2^k
                mTemp <- expand.grid(lapply(seq_along(Estim.Opt$ls.constrained.par),
                                            function(ik){return(c(0L, 1L))}))
                
                ##  Repeat the releveant columns
                mTemp <- t(Reduce(cbind, lapply(seq_along(Estim.Opt$ls.constrained.par),
                                                function(ik){
                                                    iTemp <- length(Estim.Opt$ls.constrained.par[[ik]])
                                                    mTemp <- matrix(rep(m.constraints[, ik],
                                                                        times = iTemp), ncol = iTemp)
                                                    return(mTemp)
                                                })))
                
                rownames(mTemp) <- unlist(Estim.Opt$ls.constrained.par)
                
                ##  Sort based on rownames to match parameter vector
                vSort <- rep(NA, nrow(mTemp))
                for(k in seq_len(nrow(mTemp))){
                    vSort[k] <- which(row.names(mTemp)[k] == strT)
                }

                ##  NVAR x 2^k
                mTemp <- mTemp[order(vSort), ]
            } else {
                mT <- Estim.Opt$m.constraints
                if(nrow(mT) < length(strT)) stop("Check the dimensions of the supplied matrix of constraints! \n")
            }
            
            ##  If we are working with a mixed logit -- return as a list
            if(length(Estim.Opt$ls.rand.par) > 0L){
                lsTemp <- as.list(as.data.frame(mTemp))
                strTemp <- unlist(Estim.Opt$ls.rand.par)
                lsTemp <- lapply(lsTemp, function(vX){
                    names(vX) <- strTemp
                    return(vX)
                })
                assign("ls.constraints", lsTemp, envir = .GlobalEnv)
            } else {
                assign("m.constraints", mTemp, envir = .GlobalEnv)
            }
        }
    }
    
    ############################################################################
    ##  Create a list containing the draws - IND*DRAWS x NVAR
    ############################################################################
    if(Estim.Opt$b.make.draws){
        ##  Create the matrix of draws and split to a list with elements equal
        ##  to the number of cores
        ##  IND*DRAWS x NVAR
        m.draws <- fn.generate.draws(Estim.Opt)
        ls.tmp <- vector(mode = "list", Estim.Opt$i.cores)
        ls.tmp <- lapply(seq_along(ls.id), function(i.x){
            if(i.x == 1L) i.index.end <<- 0L
            i.index <- i.index.end + (length(unique(ls.id[[i.x]])) * Estim.Opt$i.draws)
            m.tmp <- m.draws[(1L + i.index.end):(i.index), , drop = FALSE]
            i.index.end <<- i.index
            return(m.tmp)
        })
        
        ##  Check that the list elements are matrices
        ls.tmp <- lapply(ls.tmp, function(m.x) {
            if(!is.matrix(m.x)) m.x <- as.matrix(m.x)
            return(m.x)
        })
        
        ##  Save the list to global environment
        ##  IND*DRAWS x NVAR
        assign("ls.draws", ls.tmp, envir = .GlobalEnv)
        
        ##  Create a list of draws indices
        ls.tmp <- lapply(ls.tmp, function(m.x){
            ls.tmp <- split(seq_len(nrow(m.x)),
                                    rep(1L:(nrow(m.x)/Estim.Opt$i.draws),
                                        each = Estim.Opt$i.draws))
        })
        
        ##  Save the list to global environment
        ##  IND*DRAWS
        assign("ls.draws.index", ls.tmp, envir = .GlobalEnv)
        
        ########################################################################
        rm(m.draws, ls.tmp)
        rm(i.index.end, envir = .GlobalEnv)
        ########################################################################
    }
    
    ############################################################################
    rm(ls.data, ls.id, m.data, ls.data.tmp)
    invisible(gc(verbose = FALSE))
    ############################################################################
    
}

