################################################################################
##  Name:         r.script.run.lc.R
##  Created:      2018.08.08
##  Last edited:  2018.08.08
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
fn.detach.all.data <- function() {
    pos.to.detach <- (1:length(search()))[substring(search(), 
                                                    first = 1,
                                                    last = 8) != "package:" &
                                              search() != ".GlobalEnv" & 
                                              search() != "Autoloads" &
                                              search() != "CheckExEnv" &
                                              search() != "tools:rstudio" &
                                              search() != "TempEnv"]
    for (i in 1:length(pos.to.detach)) {
        if (length(pos.to.detach) > 0) {
            detach(pos = pos.to.detach[1])
            pos.to.detach <- (1:length(search()))[substring(search(), 
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
fn.detach.all.data()

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
source("r.script.ll.eclc.mnl.R")

################################################################################
##  Set up a list of options that controls how the model runs
################################################################################
Estim.Opt <- list()

################################################################################
##  Define the options for storing output
##
##  'str.output.estimates'  will be used to store .txt and .rds files of model
##                          outputs
################################################################################
Estim.Opt$str.model.name <- "Equality Constrained Latent Class - ECLC"
Estim.Opt$str.output.estimates <- "output.eclc.coral"

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
Estim.Opt$str.data <- "../data/data.coral.rds"
Estim.Opt$b.complete.data <- FALSE
Estim.Opt$i.ind <- 397L
Estim.Opt$i.obs <- 4683L
Estim.Opt$i.alts <- 3L
Estim.Opt$i.tasks <- 12L

################################################################################
##  Specify estimation options
##
##  Defaults:
##  d.tol       -   1e-06
##  d.gradtol   -   1e-06
##  d.reltol    -   1e-06
##  d.steptol   -   1e-10
##
##  str.method  -   "BFGS", "BFGSR", "NM", "CG", "NR", "BHHH", "SANN"
##
##  'b.final.hessian' is a boolean. If TRUE the final hessian from maxLik is used
##  and if FALSE the final hessian is calculated using numDeriv::hessian()
##
##  See ?numDeriv and ?maxLik
################################################################################
Estim.Opt$d.tol         <- 1e-08  
Estim.Opt$d.gradtol     <- 1e-08   
Estim.Opt$d.reltol      <- 1e-08   
Estim.Opt$d.steptol     <- 1e-10   
Estim.Opt$i.iterlim     <- 500L   
Estim.Opt$i.print.level <- 3L   
Estim.Opt$str.method    <- "BFGS" 
Estim.Opt$b.final.hessian <- FALSE

################################################################################
##  Specify estimation options for using multiple cores/parallel
##  
##  'str.output.bug'        Redirects output from workers. Useful for debugging.
##  'b.print.worker.info'   If TRUE prints information about objects loaded to 
##                          the workers
##  'i.cores'               Indicates the number of cores to use
##
##  See ?parallel
################################################################################
Estim.Opt$b.parallel <- FALSE
Estim.Opt$str.output.debug <- "" 
Estim.Opt$b.print.worker.info <- TRUE 
Estim.Opt$i.cores <- 3L
Estim.Opt$str.packages <- c("maxLik", "numDeriv", "matrixStats", "msm", "pryr",
                            "MASS")

################################################################################
##  Specify information about standard errors and scaling of utility
##
##  'b.robust.vcov'             If TRUE calculates a robust vcov using a 
##                              sandwich estimator
##  'b.adjusted.robust.vcov'    If TRUE applies and adjustment weight to the
##                              robust vcov. 
##  'b.rescale.utility'         If TRUE rescales utilitiy before taking the
##                              exponent by subtracting the middle value
################################################################################
Estim.Opt$b.robust.vcov           <- TRUE
Estim.Opt$b.adjusted.robust.vcov  <- TRUE
Estim.Opt$b.rescale.utility <- FALSE 

################################################################################
##  Specify strings including the names of the variables in your data
##
##  'str.id'        String indicating the individual id variable. This should
##                  be continuous starting at 1
##  'str.ct'        String indicating the choice task variable
##  'str.alt'       String indicating the alternative variable
##  'str.choice'    String indicating the choice variable
##  'str.cost'      String indicating the cost variable
################################################################################
Estim.Opt$str.id <- "id"
Estim.Opt$str.ct <- "ct"
Estim.Opt$str.alt <- "alt"
Estim.Opt$str.choice <- "choice"
Estim.Opt$str.cost <-  "cost"

################################################################################
##  Specify the fixed part of utility
##
##  'str.fixed.par' Character string containing the variables with non-random
##                  parameters
##  
##  If no variables have non-random parameters - leave empty: c()
################################################################################
Estim.Opt$str.fixed.par <- c("cost", "small", "large", "oil", "fish", "hab")

################################################################################
##  Estimate in willingness-to-pay space
##
##  'b.wtp.space'   If TRUE - estimate the model in WTP space, however, this 
##                  works poorly with the MNL/LC models and is sensitive to
##                  starting values.
################################################################################
Estim.Opt$b.wtp.space <- FALSE

################################################################################
##  Specify information about the draws
##
##  'b.make.draws'      If TRUE - the data set-up function will generate a matrix
##                      of draws to be used in simulation.
##  'b.correlation'     If TRUE - calculate the off-diagonal elements of the 
##                      lower Cholesky matrix. Correlation is limited to
##                      normal based distributions.
##  'str.draws.type'    String indicating the type of draws to use.
##                      Ex. 'sobol', 'halton', 'mlhs' and 'pseudo'
##  'i.drop'            Indicates the number of elements to drop from the
##                      from the sequence and is only valid for 'halton' draws.
##  'b.scramble'        If TRUE - then the sequence is scrambled. This is only
##                      valid for 'halton' and 'sobol' draws.
##  'i.scrambling.type' Indicates the type of scrambling to use for the 
##                      'sobol' draws.
##
##  See ?randtoolbox
################################################################################
Estim.Opt$b.make.draws <- FALSE
Estim.Opt$b.correlation <- FALSE
Estim.Opt$str.draws.type <- "sobol" 
Estim.Opt$i.draws <- 50L
Estim.Opt$i.drop <- 0L
Estim.Opt$b.scramble <- TRUE
Estim.Opt$i.scrambling.type <- 3L 

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
Estim.Opt$ls.rand.par <- list()

################################################################################
##  Specify a list of variables where the mean of the parameter distribution is
##  interacted with socio-demographic variables.
##
##  If no parameter distributions have interactions - leave empty: list()
##
##  'b.class.specific'  Only used for LC-MIXL. If TRUE - class specific
##                      interactoins. If FALSE - assumed to be the same across
##                      classes.
################################################################################
Estim.Opt$ls.het.par <-  list()
Estim.Opt$b.class.specific <- TRUE

################################################################################
##  Specify whether we are estimating relative scale parameters
##
##  'b.relative.scale'  If TRUE - then the model will estimate a relative
##                      scale parameter
##  'str.scale.par'     Character string containing the indicators for each
##                      group. Include all groups since the last group is
##                      automatically set to unity for identification
################################################################################
Estim.Opt$b.relative.scale <- FALSE
Estim.Opt$str.scale.par <- c()

################################################################################
##  Latent Class Model
##
##  'b.latent.class'    If TRUE - estimate a Latent Class Model 
##  'i.classes'         Indicates the number of classes (Not used for ECLC)
##  'str.class.par'     Character string containing the variables entering the
##                      class probability function. Must at the very least
##                      contain 'const', which is a generic constant added to
##                      the data if missing.
################################################################################
Estim.Opt$b.latent.class <- TRUE
Estim.Opt$i.classes <- 2L
Estim.Opt$str.class.par <- c("const")

################################################################################
##  Specify equality constraints (to use in AN-A models)
##
##  'b.equality.constrained'    If TRUE - impose an equality constraint where
##                              parameters/distributions are restricted to be
##                              equal or zero (i.e. a typical AN-A model)
##  'b.discrete.mixture           If TRUE - calculate the class probabilities as
##                              as a mixture of the probabilities of 1/0 for 
##                              each attribute. If FALSE - use MNL probs.
##  'm.constraints'             A user supplied matrix of constraints containing
##                              0 and 1. Each row an attribute and each column
##                              a class constraint. Relevant rows are repeated
##                              based on the specification in ls.constrained.par.
##                              Use same order as in 'str.fixed.par'. If a matrix
##                              is supplied then we cannot calculate probabilities
##                              as a discrete mixture. 
##                              If NULL - estimate the full 2^k. 
##  'ls.constrained.par'        List of strings containing the attributes that
##                              are constrained. In case of dummy coding and you
##                              do not consider only some levels to be constrained,
##                              you specify as a list which variables are part of
##                              the same attribute.
##                              The names of the attributes must be in the same
##                              order as ls.rand.par or fixed.par
################################################################################
Estim.Opt$b.equality.constrained <- TRUE
Estim.Opt$b.discrete.mixture <- TRUE
Estim.Opt$m.constraints <- NULL
Estim.Opt$ls.constrained.par <- list(size = c("small", "large"),
                                     industry = c("oil", "gas"),
                                     hab = c("hab"),
                                     cost = c("cost"))

################################################################################
##  Specify information about starting values
##  
##  'b.search.starting.values'  If TRUE - there will be an initial search for
##                              starting values prior to estimating the model.
##  'i.nr.of.starting.models'   Indicates the number of vectors of starting
##                              values to consider.
##  'i.nr.of.models'            Indicates the total number of models to run to
##                              completion. This is useful when checking for
##                              local optima. 
##  'd.multiplier'              is a multiplier included in the formula for 
##                              generating starting values. A larger multiplier
##                              results in larger spread of starting vectors.
##  'i.seed'                    is the seed used for the RNG. If you are using
##                              'mlhs' take note of the seed since the RNG is
##                              used to generate the sequence.
################################################################################
Estim.Opt$b.search.starting.values <- TRUE
Estim.Opt$i.nr.of.starting.models <- 1000L
Estim.Opt$i.nr.of.models <- 1L
Estim.Opt$d.multiplier <- 1.5
Estim.Opt$i.seed <- 57888385L

##  Vector of starting values the length of the number of parameters
v.param <- c(runif(6), rep(0.5, 4))

################################################################################
### Start running the model
################################################################################
fn.run.model(Estim.Opt)



