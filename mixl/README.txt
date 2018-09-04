################################################################################
##  Name:         README.txt
##  Created:      2018.08.02
##  Last edited:  2018.09.04
##  Author:       Erlend Dancke Sandorf
################################################################################

################################################################################
##  Mixed logit model - MIXL
################################################################################

The codes for the mixed logit model allows for the estimation of the model in both preference and willingness-to-pay space with correlated or uncorrelated distributions. You can specify interactions between socio-demographic variables and the means of the parameter distributions. You can specify a wide variety of draws and distributions to be used. 

################################################################################
##  Notes on use and estimation
################################################################################

- Keep the folder structure the same and open the .RProj in this folder. This will set the current folder to the working directory. Scripts used in the estimation of the model is sourced relative to this folder. 

- Currently the data needs to be supplied as a .rds file in the "data"-folder. The data needs to be in long format. The individual id variable should be numeric and start at 1 without gaps, i.e. 1:N. The alternative and choice task variables should also start indexing at 1, i.e 1:J and 1:T. See "data.coral.rds" for an example. If you want to use data stored in other formats you will need to change the corresponding line in the fn.setup.data() in the "r.script.methods.data.R"" -file. 

- You specify your model by changing the options in "r.script.run.mixl.R". This is the only file you need to change to run your model. The last line of this file: "fn.run.model()", will initiate the run sequence. 

- The procedure for generating scrambled Halton draws is very slow when a large number of draws are considered. The procedure for scrambled Sobol sequences does not work for very large scale problems (the scrambling procedure in the "randtoolbox" package will return values greater than or equal to unity). 

- A completed model will store the model object as a .rds file and the output as a .txt file for easy inspection and use in post-estimation.

- The model may be sensitive to starting values and prudence dictates to do robustness checks with respect to the vector of starting values used. 

################################################################################
##  References
################################################################################
