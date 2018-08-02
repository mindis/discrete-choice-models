################################################################################
##  Name:         r.script.methods.draws.R
##  Created:      2018.08.02
##  Last edited:  2018.08.02
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Function for generating MLHS draws
################################################################################
fn.shuffle <- function(v.inv){
    m.out <- v.inv[rank(runif(length(v.inv)))]
    return(m.out)
}

fn.generate.mlhs <- function(i.n, i.d, i.i){
    v.tmp <- seq(0, i.n - 1) / i.n
    m.out  <- matrix(0, i.n * i.i, i.d)
    
    i.j <- 1
    i.k <- 1
    
    while(i.j < i.i + 1) {
        i.k <- 1
        while(i.k < i.d + 1) {
            m.out[(1 + i.n * (i.j - 1)):(i.n * i.j), i.k] <- fn.shuffle(v.tmp + runif(1)/i.n)
            i.k <- i.k + 1
        }
        i.j <- i.j + 1
    }
    return(m.out)
}

################################################################################
##  Function for generating scrambled Halton sequences
################################################################################
fn.start <- function(i.d, i.max.digit){
    m.count <- matrix(0, i.d, i.max.digit)
    m.count[, 1] <- rep(1, i.d)
    v.digit <- rep(1, i.d)
    
    return(list(m.count, v.digit))
}

##  Function for expanding the sequence of integers (Equation 1 in Bhat(2003))
fn.digitize <- function(i.d, v.primes, m.count, v.digit){
    i.m <- 1L
    i.x <- NULL
    
    while(i.m <= i.d){
        i.l <- 1L
        i.r <- 1L
        while(i.r == 1L){
            i.x <- m.count[i.m, i.l] != (v.primes[i.m] - 1L)
            i.r <- i.r - i.x
            m.count[i.m, i.l] <- (m.count[i.m, i.l] + 1L) * (i.x == 1L)
            v.digit[i.m] <- ((i.l - 1L) == v.digit[i.m]) + v.digit[i.m]
            i.l <- i.l + 1L
        }
        i.m <- i.m + 1L
    }
    return(list(m.count, v.digit))
}

##  Function for computing the radical inverse (Equation 2 in Bhat(2003))
fn.rad.inverse <- function(i.d, v.primes, m.count, v.digit, m.perms){
    i.m <- 1L
    v.g <- matrix(0L, 1L, i.d)
    
    while(i.m <= i.d){
        i.l <- 1L
        d.p <- v.primes[i.m]
        while(i.l <= v.digit[i.m]){
            v.g[i.m] <- (m.perms[i.m, (m.count[i.m, i.l] + 1L)] / d.p) + v.g[i.m]
            d.p <- d.p * v.primes[i.m]
            i.l <- i.l + 1L
        }
        i.m <- i.m + 1L
    }
    return(v.g)
}

##  Code for generating the Halton sequence
fn.generate.halton <- function(i.n, i.d){
    if(i.d > 16L) stop("Cannot scramble Halton sequences beyond dimension 16.")
    
    i.max.digit <- 50L
    v.primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53)
    m.perms <- readRDS("../methods/m.xbrat.rds")
    
    m.h <- matrix(0, i.n, i.d)
    v.primes <- v.primes[1:i.d]
    
    ##  Get values for the first number in the halton sequence
    ls.tmp <- fn.start(i.d, i.max.digit)
    m.count <- ls.tmp[[1]]
    v.digit <- ls.tmp[[2]]
    
    m.h[1L, ] <- fn.rad.inverse(i.d, v.primes, m.count, v.digit, m.perms)
    
    ##  For the rest of the numbers
    i.j <- 2L
    while(i.j <= i.n){
        ls.tmp <- fn.digitize(i.d, v.primes, m.count, v.digit)
        m.count <- ls.tmp[[1]]
        v.digit <- ls.tmp[[2]]
        
        m.h[i.j, ] <- fn.rad.inverse(i.d, v.primes, m.count, v.digit, m.perms)
        i.j <- i.j + 1L
    }
    return(m.h)
}

################################################################################
##  Function for generating random draws
################################################################################
fn.generate.draws <- function(Estim.Opt){
    ##  Check that the specified type of draws are reckognized by the code
    if(!Estim.Opt$str.draws.type %in% c("pseudo", "halton", "sobol", "mlhs")){
        stop("Could not recognize type of draws. Please use one of the following: \n
             'pseudo', 'halton', 'sobol' or 'mlhs'.\n")
    }
    
    ##  Check that randtoolbox is loaded
    if(!("randtoolbox" %in% .packages())) require(randtoolbox)
    
    ##  Set the number of random variables
    i.K <- length(Estim.Opt$ls.rand.par)
    
    ##  Pseudo-random draws
    if(Estim.Opt$str.draws.type == "pseudo"){
        m.draws <- matrix(runif(Estim.Opt$i.ind * Estim.Opt$i.draws * i.K),
                          nrows = (Estim.Opt$i.ind * Estim.Opt$i.draws),
                          ncol = i.K)
    }
    
    ##  Halton draws
    if(Estim.Opt$str.draws.type == "halton"){
        i.tmp <- (Estim.Opt$i.ind * Estim.Opt$i.draws) + Estim.Opt$i.drop
        if(Estim.Opt$b.scramble){
            m.draws <- fn.generate.halton(i.tmp, i.K)
        } else {
            m.draws <- halton(i.tmp, dim = i.K)
        }
    }
    
    ##  Sobol draws
    if(Estim.Opt$str.draws.type == "sobol"){
        i.tmp <- Estim.Opt$i.ind * Estim.Opt$i.draws
        if(Estim.Opt$b.scramble){
            m.draws <- sobol(i.tmp, dim = i.K,
                             scrambling = Estim.Opt$i.scrambling.type)
        } else {
            m.draws <- sobol(i.tmp, dim = i.K)
        }
    }
    
    ##  MLHS draws
    if(Estim.Opt$str.draws.type == "mlhs"){
        m.draws <- fn.generate.mlhs(Estim.Opt$i.ind, i.K, Estim.Opt$i.draws)
    }
    
    ##  Return the matrix of draws
    ## IND*DRAWS x NVAR
    m.draws <- as.matrix(m.draws)
    if((Estim.Opt$str.draws.type %in% c("halton") && Estim.Opt$i.drop > 0L)){
        return(as.matrix(qnorm(m.draws[-(1:Estim.Opt$i.drop), ])))
    } else {
        return(as.matrix(qnorm(m.draws)))
    }
}



