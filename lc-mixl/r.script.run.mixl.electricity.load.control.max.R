################################################################################
################################################################################
### Name:         r.script.run.mixl.R
### Created:      2018.05.29
### Last edited:  2018.07.19
### Author:       Erlend Dancke Sandorf
### Contributors: N/A
################################################################################
################################################################################

################################################################################
### Clean memory and load functions and packages
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

##  matrixStats, parallel, maxLik, numDeriv, randtoolbox, sandwich

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
source("r.script.ll.lc.mixl.R")

################################################################################
### Set up a list of options that controls how the model runs
################################################################################
Estim.Opt <- list()

################################################################################
##  Define the options for storing output
################################################################################
Estim.Opt$str.model.name <- "Mixed Logit Model -- MIXL"
Estim.Opt$str.output.estimates <- "output.mixl.electricity.load.control.max"

################################################################################
##  Specify information about the data
##  NOTE: complete data function uses the information on i.ind, i.alt and i.task
################################################################################
Estim.Opt$str.data <- "../data/data.electricity.load.control.max.rds"
Estim.Opt$b.complete.data <- TRUE #  If FALSE not all individuals have answered all choice tasks
Estim.Opt$i.ind <- 2014L
Estim.Opt$i.obs <- 16112L
Estim.Opt$i.alts <- 3L
Estim.Opt$i.tasks <- 8L

################################################################################
##  Specify estimation options
################################################################################
Estim.Opt$d.tol         <- 1e-08   #  Tolerance level for convergence default: 1e-06
Estim.Opt$d.gradtol     <- 1e-08   #  Tolerance level for the gradient default: 1e-06
Estim.Opt$d.reltol      <- 1e-08   #  Relative tolerance for convergance default: 1e-06
Estim.Opt$d.steptol     <- 1e-08   #  Step tolerance default: 1e-10
Estim.Opt$i.iterlim     <- 500L   #  Iteration limit (zero for convergence at starting values)
Estim.Opt$i.print.level <- 3L   
Estim.Opt$str.method    <- "BFGS" # "BFGS", "BFGSR" "NM", "CG", "NR", "BHHH",
                                    # "SANN"
Estim.Opt$b.final.hessian <- FALSE ##   If FALSE don't calculate finalHessian and use numDeriv Hessian (takes longer, but is more precise) Not set up yet


################################################################################
##  Specify estimation options for parallel computing
################################################################################
Estim.Opt$b.parallel <- TRUE
Estim.Opt$str.output.debug <- "" ## Leave blank to avoid redirecting output
Estim.Opt$b.print.worker.info <- TRUE ## Currently not set up to work 
Estim.Opt$i.cores <- 3L ## detectCores() - 1
Estim.Opt$str.packages <- c("maxLik", "numDeriv", "matrixStats", "msm", "pryr",
                            "MASS")

################################################################################
##  Specify variables
################################################################################
Estim.Opt$str.id <- "id"
Estim.Opt$str.ct <- "ct.order"
Estim.Opt$str.alt <- "alt"
Estim.Opt$str.choice <- "choice"
Estim.Opt$str.cost <-  "comp.s"

################################################################################
##  Estimate in willingness-to-pay space
################################################################################
Estim.Opt$b.wtp.space <- FALSE

################################################################################
##  Specify information about the draws
################################################################################
Estim.Opt$b.make.draws <- TRUE
Estim.Opt$b.correlation <- FALSE
Estim.Opt$str.draws.type <- "sobol" # "pseudo", "halton", "sobol", "mlhs"
Estim.Opt$i.draws <- 50L
Estim.Opt$i.drop <- 0L
Estim.Opt$b.scramble <- TRUE
Estim.Opt$i.scrambling.type <- 3L ##    Only affects Sobol draws

################################################################################
##  Specify information about standard errors and scaling of utility
################################################################################
Estim.Opt$b.robust.vcov           <- TRUE
Estim.Opt$b.adjusted.robust.vcov  <- TRUE
Estim.Opt$b.rescale.utility <- TRUE # rescaling based on subraction of middle value

################################################################################
##  Specify fixed parameters
################################################################################
Estim.Opt$str.fixed.par <- c()

################################################################################
##  Random parameters
##  "n"       - Normal                                                                  
##  "ln"      - Log-normal                                                              
##  "-ln"     - Log-normal with sign change                                             
##  "u"       - Uniform -1, 1                                                           
##  "-cu"     - Constrained negative uniform                                            
##  "t"       - Triangular                                                              
##  "ct"      - Constrained triangular
##  "sb"      - Johnsen SB
##
##  Ex. <- list(att1 = "n", att2 = "-ln", att3 = "t")
##  If no interactions, leave empty: <- list()
################################################################################
Estim.Opt$ls.rand.par <- list(comp.s = "ln", max.3500 = "n", max.2000 = "n",
                              flex.use = "n", duration.90 = "n", duration.180 = "n",
                              days.10 = "n", days.20 = "n", opt.out = "n")

################################################################################
##  Specify a list of interaction variables
##  Ex. <- list(att1 = c("dem1"), att2 = c("dem1", "dem2"))
##  If no interactions, leave empty: <- list()
################################################################################
Estim.Opt$ls.het.par <-  list()

################################################################################
##  Specify whether we are estimating relative scale parameters
################################################################################
Estim.Opt$b.relative.scale <- TRUE
Estim.Opt$str.scale.par <- c("green.treatment", "brown.treatment")

################################################################################
##  Latent class model
################################################################################
Estim.Opt$b.latent.class <- TRUE
Estim.Opt$i.classes <- 2L
Estim.Opt$str.class.par <- c("const")

################################################################################
##  Specify information about starting values
################################################################################
Estim.Opt$b.search.starting.values <- TRUE
Estim.Opt$i.nr.of.starting.models <- 100L
Estim.Opt$i.nr.of.models <- 1L
Estim.Opt$d.multiplier <- 1.5
Estim.Opt$i.seed <- 57888385L

##  Vector of starting values
v.param <- c(rep(0, 18), rep(0.1, 18), 1, 0.5)

################################################################################
### Start running the model
################################################################################

fn.run.model(Estim.Opt)



