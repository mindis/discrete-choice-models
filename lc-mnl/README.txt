################################################################################
##  Name:         README.txt
##  Created:      2018.08.02
##  Last edited:  2018.09.05
##  Author:       Erlend Dancke Sandorf
################################################################################

################################################################################
##  Latent class model - LC - MNL
################################################################################

This is the code for estimating the standard latent class model. The model allows to specify the number of classes and include socio-demographic variables in the class probability function.

################################################################################
##  Notes on use and estimation
################################################################################

- Keep the folder structure the same and open the .RProj in this folder. This will set the current folder to the working directory. Scripts used in the estimation of the model is sourced relative to this folder. 

- Currently the data needs to be supplied as a .rds file in the "data"-folder. The data needs to be in long format. The individual id variable should be numeric and start at 1 without gaps, i.e. 1:N. The alternative and choice task variables should also start indexing at 1, i.e 1:J and 1:T. See "data.coral.rds" for an example. If you want to use data stored in other formats you will need to change the corresponding line in the fn.setup.data() in the "r.script.methods.data.R"" -file.

- You specify your model by changing the options in "r.script.run.lc.mnl.R". This is the only file you need to change to run your model. The last line of this file: "fn.run.model()", will initiate the run sequence. 

- A completed model will store the model object as a .rds file and the output as a .txt file for easy inspection and use in post-estimation.

- The model is sensitive to starting values and prudence dictates to do robustness checks with respect to the vector of starting values used. 

################################################################################
##  References
################################################################################
