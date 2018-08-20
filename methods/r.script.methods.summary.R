################################################################################
##  Name:         r.script.methods.summary.R
##  Created:      2018.08.02
##  Last edited:  2018.08.02
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Create a function to calculate and print summary statistics
################################################################################
fn.model.summary <- function(ls.X){
    ##  Set this option to avoid sceintific notation in printing
    options(scipen = 999)
    
    ############################################################################
    ##  Print general convergence information
    ############################################################################
    cat("###################################################################\n")
    cat("Model name: ", Estim.Opt$str.model.name, "\n")
    cat("Model message: ", ls.X$message, "\n")
    cat("Estimating model ", ls.X$i.model.number, "out of ",
        Estim.Opt$i.nr.of.models, ". (", ls.X$i.count.failed,
        "model(s) have failed).\n\n")
    
    cat("###################################################################\n")
    cat("Starting values (NOTE: Padded with zeros)\n")
    i.tmp <- length(ls.X$v.starting.values)
    v.tmp <- c(ls.X$v.starting.values, rep(0, ((5L * ceiling(i.tmp / 5L)) - i.tmp)))
    print(matrix(v.tmp, ncol = 5L, byrow = TRUE))
    cat("\n")
    
    cat("###################################################################\n")
    cat("Final values (NOTE: Padded with zeros)\n")
    i.tmp <- length(ls.X$estimate)
    v.tmp <- c(ls.X$estimate, rep(0, ((5L * ceiling(i.tmp / 5L)) - i.tmp)))
    print(matrix(v.tmp, ncol = 5L, byrow = TRUE))
    cat("\n")
    
    cat("###################################################################\n")
    cat("Gradient at convergence (NOTE: Padded with zeros)\n")
    i.tmp <- length(ls.X$gradient)
    v.tmp <- c(ls.X$gradient, rep(0, ((5L * ceiling(i.tmp / 5L)) - i.tmp)))
    print(matrix(v.tmp, ncol = 5L, byrow = TRUE))
    cat("\n")
    
    ############################################################################
    ##  Create a function for printing some summary statistics
    ############################################################################
    fn.summary.stats <- function(d.ll.final, d.ll.zero, i.k){
        i.tmp <- Estim.Opt$i.obs
        v.tmp <- c(d.ll.final,
                   d.ll.zero,
                   i.tmp,
                   i.k,
                   (1L - ((d.ll.final - i.k) / (d.ll.zero))), 
                   ((-2L * d.ll.final) + (2L * i.k) ), 
                   ((-2L * d.ll.final) + (3L * i.k)), 
                   ((-2L * d.ll.final) + (i.k * (log(i.tmp) + 1L))), 
                   ((-2L * d.ll.final) + (i.k * (log((i.tmp + 2L) / 24L) + 1L ))), 
                   ((-2L * d.ll.final) + (2L * i.k)  +
                        (((2L * (i.k + 1L)) * (i.k + 2L))/(i.tmp - i.k - 2L))),
                   ((-2L * d.ll.final) + (i.k * log(i.tmp))), 
                   ((-2L * d.ll.final) + (i.k * (log((i.tmp + 2L) / 24L)))), 
                   ((-2L * d.ll.final) + (i.k * (log(i.tmp) - log(2L * pi)))), 
                   ((-2L * d.ll.final) + (2L * (i.k * (log(log(i.tmp)))))))
        m.tmp <- matrix(v.tmp, nrow = 14L, ncol = 1L,
                        dimnames = list(c("LL", "LL(0)", "N", "K", "Adj. Rho^2",
                                          "AIC", "AIC3", "CAIC", "CAIC*",
                                          "HT-AIC/AICc", "BIC", "BIC*", "DBIC",
                                          "HQIC"), c("Value")))
        
        ##  Return the matrix
        return(m.tmp)
    }
    
    cat("###################################################################")
    cat("\n")
    cat("Model details: \n")
    m.stat <- fn.summary.stats(ls.X$maximum, ls.X$d.ll.zero, length(ls.X$estimate))
    print(m.stat)
    cat("\n")
    
    ############################################################################
    ##  Calculate the variance covariance matrix and construct table of outputs
    ############################################################################
    ##  Check that we don't have any missing values in the Hessian matrix
    if(!any(is.nan(ls.X$hessian))){
        m.vcov <- MASS::ginv(-ls.X$hessian)
        i.k <- length(ls.X$estimate)
        
        ##  Check whether we are calculating robust standard errors
        if(Estim.Opt$b.robust.vcov){
            m.bread <- m.vcov * Estim.Opt$i.ind
            m.bread[is.na(m.bread)] <- 0L
            
            m.meat <- crossprod(ls.X$gradientObs) / Estim.Opt$i.ind
            if(Estim.Opt$b.adjusted.robust.vcov){
                m.meat <- m.meat * (Estim.Opt$i.ind / (Estim.Opt$i.ind - i.k))
            }
            m.meat[is.na(m.meat)] <- 0L
            
            m.vcov <- (m.bread %*% m.meat %*% m.bread) / Estim.Opt$i.ind
        }
        
        ##  Matrix of outputs
        m.out <- matrix(NA, nrow = length(ls.X$estimate), ncol = 4L)
        m.out[, 1L] <- ls.X$estimate
        m.out[, 2L] <- sqrt(diag(m.vcov))
        m.out[, 3L] <- m.out[, 1L] / m.out[, 2L]
        m.out[, 4L] <- 2L * pt(-abs(m.out[, 3L]), df = Estim.Opt$i.obs)
        m.out <- round(m.out, 5)
        colnames(m.out) <- c("Est.", "Std. Err.", "T-stat", "P-value")
        rownames(m.out) <- names(ls.X$estimate)
        
        cat("###################################################################")
        cat("\n")
        cat("Estimated parameters: \n")
        print(m.out)
        cat("\n")
        
        ########################################################################
        ##  Print extra information if we have a MIXL with correlations
        ########################################################################
        if(Estim.Opt$b.correlation == TRUE){
            i.F <- length(Estim.Opt$str.fixed.par)
            i.R <- length(Estim.Opt$ls.rand.par)
            str.names <- names(Estim.Opt$ls.rand.par)
            i.start <- (1L + i.F + i.R)
            i.end <- (i.F + ((2L * i.R) + (i.R * (i.R - 1L) / 2L)))
            
            ##   Set up the lower Cholesky matrix
            m.L <- matrix(0, i.R, i.R)
            m.L[lower.tri(m.L, diag = TRUE)] <- ls.X$estimate[i.start:i.end]
            
            ##   Create the correlation matrix
            m.vcov.param <- m.L %*% t(m.L)
            m.D <- solve(diag(sqrt(diag(m.vcov.param))))
            m.rho <- m.D %*% m.vcov.param %*% m.D
            
            ##   Combine the lower Cholesky and upper triangular correlation matrices
            m.tmp <- matrix(0, i.R, i.R)
            m.tmp[upper.tri(m.rho, diag = FALSE)] <- m.rho[lower.tri(m.rho, diag = FALSE)]
            m.chol.corr <- m.L + m.tmp
            colnames(m.chol.corr) <- str.names
            rownames(m.chol.corr) <- str.names
            
            ##   Subset the vcov to match the lower Cholesky matrix
            m.lambda <- m.vcov[i.start:i.end, i.start:i.end]
            
            ##   Calculate the SDs of the random parameters
            m.sd <- fn.sd(m.L, m.lambda)
            colnames(m.sd) <- c("Est.", "Std. Err.", "T-stat", "P-value")
            rownames(m.sd) <- str.names
            
            ##   Calculate the standard errors of the correlations between the parameters
            m.corr <- fn.corr(m.L, m.lambda)
            colnames(m.corr) <-c("Est.", "Std. Err.", "T-stat", "P-value")
            str.corr.names <- NULL
            
            for(i in 1:(length(Estim.Opt$ls.rand.par) - 1)){
                str.corr.names <- c(str.corr.names, paste(str.names[i],
                                          str.names[(i + 1):i.R], sep = "."))
            }
            
            rownames(m.corr) <- str.corr.names
            
            cat("###################################################################")
            cat("\n")
            cat("Standard deviations of the random parameters \n")
            print(m.sd)
            cat("\n")
            
            cat("###################################################################")
            cat("\n")
            cat("Lower triangular Cholesky matrix and upper triangular correlation matrix \n")
            print(m.chol.corr)
            cat("\n")
            
            cat("###################################################################")
            cat("\n")
            cat("Correlations between the random parameters \n")
            print(m.corr)
            cat("\n")
        }
        
        ########################################################################
        ##  Print the matrix of restrictions if we are estimating the ECLC models
        ########################################################################
        if(Estim.Opt$b.equality.constrained){
            str.col.names <- paste0("class.", seq_len(ncol(m.delta.expanded)))
            m.tmp <- m.delta.expanded
            colnames(m.tmp) <- str.col.names
            cat("###################################################################")
            cat("\n")
            cat("The class restrictions: \n")
            print(m.tmp)
            cat("\n")
        }
    } else {
        cat("###################################################################")
        cat("\n")
        cat("The Hessian contains NA. Moving to next model..")
        cat("\n")
    }
    
    
    tmp = fn.time(ls.X$time)
    cat("###################################################################")
    cat("\n")
    cat("Model was estimated on", ifelse(Estim.Opt$b.parallel == TRUE, Estim.Opt$i.cores, "1"), "core(s).\n")
    cat("Model completion time:", tmp[[1]], "days", tmp[[2]],"hours", tmp[[3]], "minutes and", tmp[[4]], "seconds \n")
    if(Estim.Opt$b.make.draws){
        cat("Estimated", length(Estim.Opt$ls.rand.par), "random variables using", Estim.Opt$i.draws, ifelse(Estim.Opt$b.scramble == TRUE, "scrambled", ""), Estim.Opt$str.draws.type, "draws.\n")
    }
    cat("The model is estimated in", ifelse(Estim.Opt$b.wtp.space == TRUE, "WTP space.\n", "preference space.\n"))
    cat("The results are shown with", ifelse(Estim.Opt$b.adjusted.robust.vcov == TRUE & Estim.Opt$b.robust.vcov == TRUE, "adjusted", ""), ifelse(Estim.Opt$b.robust.vcov == TRUE, "robust", "normal"), "standard errors.\n")
    cat("The model used the following seed: ", ls.X$i.seed, "\n")
    cat(Estim.Opt$i.ind, "respondents made", Estim.Opt$i.obs, "choices.\n")
    cat("###################################################################")
    cat("\n")
}

################################################################################
##  Function for calculating time spent running the model
################################################################################
fn.time = function(t) {
    if(class(t) != "proc_time") stop("Input must be of class 'proc_time'")
    
    i.t <- as.numeric(t[3])
    
    i.days    <- floor( i.t / 86400)
    i.hours   <- floor((i.t - (i.days * 86400)) / 3600)
    i.minutes <- floor((i.t - (i.days * 86400) - (i.hours * 3600)) / 60)
    i.seconds <- floor( i.t - (i.days * 86400) - (i.hours * 3600) - (i.minutes * 60))
    
    return(list(i.days, i.hours, i.minutes, i.seconds))
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
            x[n, 4] = 2 * pt(-abs(y1 / y2), df = Estim.Opt$i.obs)
            
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
        x[i, 4] = 2 * pt(-abs(y1 / y2), df = Estim.Opt$i.obs)
        
        i = i + 1
        
    }
    return(x)
}




