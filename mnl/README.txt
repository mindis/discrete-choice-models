################################################################################
##  Name:         README.txt
##  Created:      2018.08.02
##  Last edited:  2018.09.05
##  Author:       Erlend Dancke Sandorf
################################################################################

################################################################################
##  Multinomial Logit Model - MNL
################################################################################

Standard multinomial logit model.

################################################################################
##  Notes on use and estimation
################################################################################

- Keep the folder structure the same and open the .RProj in this folder. This will set the current folder to the working directory. Scripts used in the estimation of the model is sourced relative to this folder. 

- Currently the data needs to be supplied as a .rds file in the "data"-folder. The data needs to be in long format. The individual id variable should be numeric and start at 1 without gaps, i.e. 1:N. The alternative and choice task variables should also start indexing at 1, i.e 1:J and 1:T. See "data.coral.rds" for an example. If you want to use data stored in other formats you will need to change the corresponding line in the fn.setup.data() in the "r.script.methods.data.R"" -file. 

- You specify your model by changing the options in "r.script.run.mnl.R". This is the only file you need to change to run your model. The last line of this file: "fn.run.model()", will initiate the run sequence. 

- A completed model will store the model object as a .rds file and the output as a .txt file for easy inspection and use in post-estimation.

################################################################################
##  References
################################################################################

