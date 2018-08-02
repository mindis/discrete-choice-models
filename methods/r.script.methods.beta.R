################################################################################
##  Name:         r.script.methods.beta.R
##  Created:      2018.08.02
##  Last edited:  2018.08.02
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for generating beta names
################################################################################
fn.make.beta.names <- function(Estim.Opt){
    ##  Create empty placeholders for the strings of parameter names
    str.fixed <- str.mean <- str.std <- str.het <- str.scale <- str.class <- NULL
    
    ##  Fixed parameters
    if(length(Estim.Opt$str.fixed.par) > 0){
        str.fixed <- paste0("beta.", Estim.Opt$str.fixed.par)
        if(Estim.Opt$b.latent.class){
            str.fixed <- paste("beta", rep(Estim.Opt$str.fixed.par,
                                           times = Estim.Opt$i.classes),
                               rep(seq_len(Estim.Opt$i.classes),
                                   each = length(Estim.Opt$str.fixed.par)),
                               sep = ".")
        }
    }
    
    ##  Random parameters
    if(length(Estim.Opt$ls.rand.par) > 0){
        str.tmp <- names(Estim.Opt$ls.rand.par)
        str.mean <- paste0("beta.", str.tmp)
        str.std <- paste0("sigma.", str.tmp)
        if(Estim.Opt$b.correlation){
            str.std <- NULL
            i.K <- length(Estim.Opt$ls.rand.par)
            for(i in 1:i.K){
                str.std <- c(str.std,
                             paste("sigma", str.tmp[i], str.tmp[i:i.K],
                                   sep = "."))
            }
        }
    }
    
    ##  Interactions with random parameters
    if(length(Estim.Opt$ls.het.par) > 0 ){
        str.tmp <- names(Estim.Opt$ls.het.par)
        str.het <- unlist(lapply(str.tmp, function(x){
            outer(x, Estim.Opt$ls.het.par[[x]], FUN = paste, sep = ".")
        }))
        str.het <- paste0("phi.", str.het)
    }
    
    ##  Relative scale parameters
    if(Estim.Opt$b.relative.scale){
        i.K <- length(Estim.Opt$str.scale.par)
        str.scale <- paste0("lambda.", Estim.Opt$str.scale.par[-c(i.K)])
    }
    
    ## Class probability function parameters
    if(Estim.Opt$b.latent.class){
        i.K <- length(Estim.Opt$str.class.par)
        i.Q <- Estim.Opt$i.classes
        str.class <- paste("theta",
                           rep(Estim.Opt$str.class.par, times = (i.Q - 1L)),
                           rep(seq_len((i.Q - 1L)), each = i.K),
                           sep = ".")
    }
    
    ##  Return the list of parameter names
    return(list(str.fixed = str.fixed,
                str.mean = str.mean,
                str.std = str.std,
                str.het = str.het,
                str.scale = str.scale,
                str.class = str.class))
}

################################################################################
##  Function for generating betas
################################################################################
fn.make.beta <- function(Estim.Opt, v.beta, v.sigma, m.draws){
    ##  Check that the specified distributions are valid
    if(any(Estim.Opt$ls.rand.par %in% c("n", "ln", "-ln",
                                        "u", "-cu", "t",
                                        "ct", "sb") == FALSE)){
        stop("Cannot recognize the specified distributions.")
    }
    
    ##  IND*DRAWS x NVAR
    colnames(m.draws) <- names(Estim.Opt$ls.rand.par)
    
    ############################################################################
    ##  Function to transform the distributions
    ############################################################################
    fn.transform.distributions <- function(i.x){
        ##  Get the position of the draws for the corresponding beta
        # i.pos <- grep(names(Estim.Opt$ls.rand.par[i.x]), names(v.beta))
        i.pos <- which(names(v.beta) == paste0("beta.",
                                               names(Estim.Opt$ls.rand.par)[i.x]))
        
        ##  Negative log-normal
        if(Estim.Opt$ls.rand.par[[i.x]] == "-ln"){
            if(Estim.Opt$b.correlation){
                v.tmp <- exp(m.beta[, i.pos]) * -1
            } else {
                v.tmp <- exp(v.beta[i.pos] + v.sigma[i.pos] * m.draws[, i.pos, drop = FALSE]) * -1
            }
        }
        
        ##  Log-normal
        if(Estim.Opt$ls.rand.par[[i.x]] == "ln"){
            if(Estim.Opt$b.correlation){
                v.tmp <- exp(m.beta[, i.pos])
            } else {
                v.tmp <- exp(v.beta[i.pos] + v.sigma[i.pos] * m.draws[, i.pos, drop = FALSE])
            }
        }
        
        ##  Normal
        if(Estim.Opt$ls.rand.par[[i.x]] == "n"){
            if(Estim.Opt$b.correlation){
                v.tmp <- m.beta[, i.pos]
            } else {
                v.tmp <- v.beta[i.pos] + v.sigma[i.pos] * m.draws[, i.pos, drop = FALSE]
            }
        }
        
        ##  Uniform -1, 1
        if(Estim.Opt$ls.rand.par[[i.x]] == "u"){
            v.tmp.draws <- 2 * pnorm(m.draws[, i.pos, drop = FALSE]) - 1
            v.tmp <- v.beta[i.pos] + v.sigma[i.pos] * v.tmp.draws
        }
        
        ##  Constrained negative uniform
        if(Estim.Opt$ls.rand.par[[i.x]] == "-cu"){
            v.tmp.draws <- 2 * pnorm(m.draws[, i.pos, drop = FALSE]) - 1
            v.tmp <- exp(v.beta[i.pos] + v.sigma[i.pos] * v.tmp.draws) * -1
        }
        
        ##  Triangular
        if(Estim.Opt$ls.rand.par[[i.x]] == "t"){
            v.tmp.draws <- pnorm(m.draws[, i.pos, drop = FALSE])
            v.b.tmp.draws <- v.tmp.draws < 0.5
            v.tmp.draws <- v.b.tmp.draws * (sqrt(2 * v.tmp.draws) - 1) +
                !v.b.tmp.draws * (1 - sqrt(2 * (1 - v.tmp.draws)))
            v.tmp <- v.beta[i.pos] + v.sigma[i.pos] * v.tmp.draws
        }
        
        ##  Constrained triangular
        if(Estim.Opt$ls.rand.par[[i.x]] == "ct"){
            v.tmp.draws <- pnorm(m.draws[, i.pos, drop = FALSE])
            v.b.tmp.draws <- v.tmp.draws < 0.5
            v.tmp.draws <- v.b.tmp.draws * (sqrt(2 * v.tmp.draws) - 1) +
                !v.b.tmp.draws * (1 - sqrt(2 * (1 - v.tmp.draws)))
            d.spread <- v.beta[i.pos] * (1 / (1 + exp(v.sigma[i.pos])))
            v.tmp <- v.beta[i.pos] + d.spread * v.tmp.draws
        }
        
        ##  Johnson SB -- Not set up to work properly
        # if(Estim.Opt$ls.rand.par[[i.x]] == "sb"){
        #     v.tmp <- v.beta[i.pos] + v.sigma[i.pos] * m.draws[, i.pos, drop = FALSE]
        #     v.tmp <- exp(v.tmp) / (1 + exp(v.tmp))
        # }
        
        ##  Return the transformed parameter
        return(v.tmp)
    }
    
    ############################################################################
    ##  Check for correlation
    ############################################################################
    if(Estim.Opt$b.correlation){
        if(any(Estim.Opt$ls.rand.par %in% c("-ln", "ln", "n") == FALSE)){
            stop("Only distributions in the family of normals can be correlated")
        }
        
        ##  Set up the lower Cholesky matrix
        i.K <- length(Estim.Opt$ls.rand.par)
        m.L <- matrix(0L, i.K, i.K)
        m.L[lower.tri(m.L, diag = TRUE)] <- v.sigma
        
        ##  Calculate the betas
        m.sigma <- tcrossprod(m.L, m.draws)
        
        ##  IND*DRAWS x NVAR
        m.beta <- v.beta + m.sigma
        m.beta <- t(m.beta)
    } else {
        ##  Transform the distributions
        ls.beta <- lapply(seq_along(Estim.Opt$ls.rand.par), function(i.x){
            v.x <- fn.transform.distributions(i.x)
        })
        
        m.beta <- Reduce(cbind, ls.beta)
    }

    ## IND*DRAWS x NVAR
    if(!is.matrix(m.beta)) m.beta <- as.matrix(m.beta)
    colnames(m.beta) <- names(Estim.Opt$ls.rand.par)
    return(m.beta)
}
