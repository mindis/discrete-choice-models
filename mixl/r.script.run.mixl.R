################################################################################
##  Name:         r.script.run.mixl.R
##  Created:      2018.08.02
##  Last edited:  2018.09.11
##  Author:       Erlend Dancke Sandorf
##  Contributors: N/A
################################################################################

################################################################################
##  Clean memory and load functions and packages
################################################################################
rm(list = ls(all = TRUE))
options(prompt = "R> ", digits = 6)

################################################################################
## Function for detaching all data from all environments
################################################################################
fnDetachAllData <- function() {
    iP <- (1:length(search()))[substring(search(), 
                                         first = 1,
                                         last = 8) != "package:" &
                                   search() != ".GlobalEnv" & 
                                   search() != "Autoloads" &
                                   search() != "CheckExEnv" &
                                   search() != "tools:rstudio" &
                                   search() != "TempEnv"]
    for (i in 1:length(iP)) {
        if (length(iP) > 0) {
            detach(pos = iP[1])
            iP <- (1:length(search()))[substring(search(), 
                                                 first = 1,
                                                 last = 8) != "package:" &
                                           search() != ".GlobalEnv" &
                                           search() != "Autoloads" &
                                           search() != "CheckExEnv" &
                                           search() != "tools:rstudio" & 
                                           search() != "TempEnv"]
        }
    }
}

##  Detach all data
fnDetachAllData()

##  Turn off file-writing if in use
if(sink.number() > 0) sink()

##   Load packages
invisible(lapply(c("maxLik", "sandwich", "numDeriv", "randtoolbox", "matrixStats",
                   "msm", "parallel", "pryr", "compiler", "MASS"),
                 require, character.only = TRUE))

################################################################################
##  Load the other functions and files
################################################################################
source("../methods/r.script.methods.beta.R")
source("../methods/r.script.methods.data.R")
source("../methods/r.script.methods.draws.R")
source("../methods/r.script.methods.miscellaneous.R")
source("../methods/r.script.methods.run.R")
source("../methods/r.script.methods.summary.R")
source("../methods/r.script.methods.parallel.R")
source("r.script.ll.mixl.R")

################################################################################
##  Set up a list of options that controls how the model runs
################################################################################
EstimOpt <- list()

################################################################################
##  Define the options for storing output
##
##  'str.output.estimates'  will be used to store .txt and .rds files of model
##                          outputs
################################################################################
EstimOpt$strModelName <- "Mixed Logit Model -- MIXL"
EstimOpt$strOutputFile <- "output.mixl.coral"

################################################################################
##  Specify information about the data
##
##  'b.complete.data' is a boolean indicating whether all respondents answered
##  all choice tasks
##
##  'i.ind'     indicates the number of individuals/respondents
##  'i.obs'     indicates the number of choice observations
##  'i.alts'    indicates the number of alternatives
##  'i.tasks'   indicates the number of choice tasks
##
##  The complete data function uses the information on i.ind, i.alt and i.task
################################################################################
EstimOpt$strData <- "../data/data.coral.rds"
EstimOpt$bCompleteData <- FALSE
EstimOpt$iN <- 397L
EstimOpt$iN_obs <- 4683L
EstimOpt$iJ <- 3L
EstimOpt$iT <- 12L

################################################################################
##  Specify estimation options
##
##  Defaults:
##  dTol       -   1e-06
##  dGradtol   -   1e-06
##  dReltol    -   1e-06
##  dSteptol   -   1e-10
##
##  strMethod  -   "BFGS", "BFGSR", "NM", "CG", "NR", "BHHH", "SANN"
##
##  'bFinalHessian' is a boolean. If TRUE the final hessian from maxLik is used
##  and if FALSE the final hessian is calculated using numDeriv::hessian()
##
##  See ?numDeriv and ?maxLik, ?maxNR, ?optim
################################################################################
EstimOpt$dTol         <- 1e-08  
EstimOpt$dGradtol     <- 1e-08   
EstimOpt$dReltol      <- 1e-08   
EstimOpt$dSteptol     <- 1e-10   
EstimOpt$iIterlim     <- 500L   
EstimOpt$iPrintLevel <- 3L   
EstimOpt$strMethod    <- "BFGS" 
EstimOpt$bFinalHessian <- FALSE

################################################################################
##  Specify estimation options for using multiple cores/parallel
##  
##  'strDebugFile'        Redirects output from workers. Useful for debugging.
##  'bPrintWorkerInfo'   If TRUE prints information about objects loaded to 
##                          the workers
##  'iCores'               Indicates the number of cores to use
##
##  See ?parallel
################################################################################
EstimOpt$bParallel <- FALSE
EstimOpt$strDebugFile <- "" 
EstimOpt$bPrintWorkerInfo <- TRUE 
EstimOpt$iCores <- 3L
EstimOpt$strPackages <- c("maxLik", "numDeriv", "matrixStats", "msm", "pryr",
                          "MASS")

################################################################################
##  Specify information about standard errors and scaling of utility
##
##  'bRobustVCOV'             If TRUE calculates a robust vcov using a 
##                            sandwich estimator
##  'bRobustVCOV'             If TRUE applies and adjustment weight to the
##                            robust vcov. 
##  'bRescaleUtility'         If TRUE rescales utilitiy before taking the
##                            exponent by subtracting the middle value
################################################################################
EstimOpt$bRobustVCOV           <- TRUE
EstimOpt$bAdjustedRobustVCOV  <- TRUE
EstimOpt$bRescaleUtility <- FALSE

################################################################################
##  Specify strings including the names of the variables in your data
##
##  'strID'        String indicating the individual id variable. This should
##                 be continuous starting at 1
##  'strCT'        String indicating the choice task variable
##  'strALT'       String indicating the alternative variable
##  'strChoice'    String indicating the choice variable
##  'strP_cost'      String indicating the cost variable
################################################################################
EstimOpt$strID <- "id"
EstimOpt$strCT <- "ct"
EstimOpt$strALT <- "alt"
EstimOpt$strChoice <- "choice"
EstimOpt$strP_cost <-  "cost"

################################################################################
##  Specify the fixed part of utility
##
##  'strP_fixed' Character string containing the variables with non-random
##                  parameters
##  
##  If no variables have non-random parameters - leave empty: c()
################################################################################
EstimOpt$strP_fixed <- c()

################################################################################
##  Estimate in willingness-to-pay space
##
##  'bWTP_space'   If TRUE - estimate the model in WTP space, however, this 
##                  works poorly with the MNL/LC models and is sensitive to
##                  starting values.
################################################################################
EstimOpt$bWTP_space <- FALSE

################################################################################
##  Specify information about the draws
##
##  'bMakeDraws'    If TRUE - the data set-up function will generate a matrix
##                  of draws to be used in simulation.
##  'bCorrelation'  If TRUE - calculate the off-diagonal elements of the 
##                  lower Cholesky matrix. Correlation is limited to
##                  normal based distributions.
##  'strDrawType'   String indicating the type of draws to use.
##                  Ex. 'sobol', 'halton', 'mlhs' and 'pseudo'
##  'iR'            Number of draws
##  'iD'            Indicates the number of elements to drop from the
##                  from the sequence and is only valid for 'halton' draws.
##  'bScramble'     If TRUE - then the sequence is scrambled. This is only
##                  valid for 'halton' and 'sobol' draws.
##  'iScrambleType' Indicates the type of scrambling to use for the 
##                  'sobol' draws.
##
##  See ?randtoolbox
################################################################################
EstimOpt$bMakeDraws <- TRUE
EstimOpt$bCorrelation <- FALSE
EstimOpt$strDrawType <- "halton" 
EstimOpt$iR <- 10L
EstimOpt$iD <- 10L
EstimOpt$bScramble <- FALSE
EstimOpt$iScrambleType <- 3L 

################################################################################
##  Specify a list of variables that have random parameters and indicate what
##  distributions they follow. Currently the codes support the following 
##  choice of distributions
##
##  "n"       - Normal                                                                  
##  "ln"      - Log-normal                                                              
##  "-ln"     - Log-normal with sign change                                             
##  "u"       - Uniform -1, 1                                                           
##  "-cu"     - Constrained negative uniform                                            
##  "t"       - Triangular                                                              
##  "ct"      - Constrained triangular
##
##  If no variables have random parameters - leave empty: list()
################################################################################
EstimOpt$lsP_rand <- list(cost = "-ln", small = "n", large = "n",
                          oil = "n", fish = "n", hab = "n")

################################################################################
##  Specify a list of variables where the mean of the parameter distribution is
##  interacted with socio-demographic variables.
##
##  If no parameter distributions have interactions - leave empty: list()
##
##  'bClassSpecific'    Only used for LC-MIXL. If TRUE - class specific
##                      interactoins. If FALSE - assumed to be the same across
##                      classes.
################################################################################
EstimOpt$lsP_het <-  list()
EstimOpt$bClassSpecific <- FALSE

################################################################################
##  Specify whether we are estimating relative scale parameters
##
##  'bRelativeScale'    If TRUE - then the model will estimate a relative
##                      scale parameter
##  'strP_scale'        Character string containing the indicators for each
##                      group. Include all groups since the last group is
##                      automatically set to unity for identification
################################################################################
EstimOpt$bRelativeScale <- FALSE
EstimOpt$strP_scale <- c("above.median", "below.median")

################################################################################
##  Latent Class Model
##
##  'bLatentClass'    If TRUE - estimate a Latent Class Model
##  'iQ'              Indicates the number of classes (Not used for ECLC)
##  'strP_class'      Character string containing the variables entering the
##                    class probability function. Must at the very least
##                    contain 'const', which is a generic constant added to
##                    the data if missing.
################################################################################
EstimOpt$bLatentClass <- FALSE
EstimOpt$iQ <- 2L
EstimOpt$strP_class <- c("const")

################################################################################
##  Specify equality constraints (to use in AN-A models)
##
##  'bEqualityConstrained'  If TRUE - impose an equality constraint where
##                          parameters/distributions are restricted to be
##                          equal or zero (i.e. a typical AN-A model)
##  'bDiscreteMixture       If TRUE - calculate the class probabilities as
##                          as a mixture of the probabilities of 1/0 for 
##                          each attribute. If FALSE - use MNL probs.
##  'mConstraints'          A user supplied matrix of constraints containing
##                          0 and 1. Each row an attribute and each column
##                          a class constraint. Relevant rows are repeated
##                          based on the specification in ls.constrained.par.
##                          Use same order as in 'str.fixed.par'. If a matrix
##                          is supplied then we cannot calculate probabilities
##                          as a discrete mixture. 
##                          If NULL - estimate the full 2^k. 
##  'lsP_constrained'       List of strings containing the attributes that
##                          are constrained. In case of dummy coding and you
##                          do not consider only some levels to be constrained,
##                          you specify as a list which variables are part of
##                          the same attribute.
##                          The names of the attributes must be in the same
##                          order as ls.rand.par or fixed.par
################################################################################
EstimOpt$bEqualityConstrained <- FALSE
EstimOpt$bDiscreteMixture <- TRUE
EstimOpt$mConstraints <- NULL 
EstimOpt$lsP_constrained <- list()

################################################################################
##  Specify information about starting values
##  
##  'bSearchStartingValues'     If TRUE - there will be an initial search for
##                              starting values prior to estimating the model.
##  'iM_search'                 Indicates the number of vectors of starting
##                              values to consider.
##  'iM_run'                    Indicates the total number of models to run to
##                              completion. This is useful when checking for
##                              local optima. 
##  'dMulti'                    is a multiplier included in the formula for 
##                              generating starting values. A larger multiplier
##                              results in larger spread of starting vectors.
##  'iSeed'                     is the seed used for the RNG. If you are using
##                              'mlhs' take note of the seed since the RNG is
##                              used to generate the sequence.
################################################################################
EstimOpt$bSearchStartingValues<- FALSE
EstimOpt$iM_search <- 100L
EstimOpt$iM_run <- 1L
EstimOpt$dMulti <- 1.5
EstimOpt$iSeed <- 57888385L

##  Vector of starting values the length of the number of parameters
vP <-  c(rep(0, 6), rep(0.1, 6))

################################################################################
### Start running the model
################################################################################
fnRunModel(EstimOpt)



