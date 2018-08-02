################################################################################
##  Name:         r.script.run.lc.R
##  Created:      2018.08.02
##  Last edited:  2018.08.02
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
source("r.script.ll.lc.R")

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
Estim.Opt$str.model.name <- "Latent Class Model - LC"
Estim.Opt$str.output.estimates <- "output.lc.demo"

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
Estim.Opt$str.data <- "../data/data.demo.rds"
Estim.Opt$b.complete.data <- TRUE
Estim.Opt$i.ind <- 200L
Estim.Opt$i.obs <- 2400L
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
Estim.Opt$str.cost <-  "x4"

################################################################################
##  Specify the fixed part of utility
##
##  'str.fixed.par' Character string containing the variables with non-random
##                  parameters
##  
##  If no variables have non-random parameters - leave empty: c()
################################################################################
Estim.Opt$str.fixed.par <- c("x4", "x1", "x2", "x3")

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
################################################################################
Estim.Opt$ls.het.par <-  list()

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
##  'i.classes'         Indicates the number of classes
##  'str.class.par'     Character string containing the variables entering the
##                      class probability function. Must at the very least
##                      contain 'const', which is a generic constant added to
##                      the data if missing.
################################################################################
Estim.Opt$b.latent.class <- TRUE
Estim.Opt$i.classes <- 3L
Estim.Opt$str.class.par <- c("const")

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
Estim.Opt$b.search.starting.values <- FALSE
Estim.Opt$i.nr.of.starting.models <- 100L
Estim.Opt$i.nr.of.models <- 1L
Estim.Opt$d.multiplier <- 1.5
Estim.Opt$i.seed <- 57888385L

##  Vector of starting values the length of the number of parameters
v.param <- c(rep(0, 12), 0.5, 0.5)

################################################################################
### Start running the model
################################################################################
fn.run.model(Estim.Opt)



