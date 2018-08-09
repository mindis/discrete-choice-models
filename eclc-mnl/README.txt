################################################################################
##  Name:         README.txt
##  Created:      2018.08.09
##  Last edited:  2018.08.09
##  Author:       Erlend Dancke Sandorf
################################################################################

################################################################################
##  Equality Constrained Latent Class Model - ECLC MNL
################################################################################

This is the equality constrained latent class model (ECLC) implemented to infer attribute non-attendance. The underlying preference structure is assumed to be a MNL model. By default the model will consider the full 2^k AN-A patterns, but the user can supply a matrix of restrictions to consider a subset. This matrix must be of the correct dimensions. Few checks are built in to ensure that this is the case and it is up to the user to do this correctly. The class probabilities can be estimated using a MNL model OR if the option 'b.discrete.mixture' is TRUE, a discrete mixture distribution will be used and the model is equivalen to the Endogenous Attribute Non-Attendance model developed by Arne Risa Hole. 

################################################################################
##  Notes on use and estimation
################################################################################

- Keep the folder structure the same and open the .RProj in this folder. This will set the current folder to the working directory. Scripts used in the estimation of the model is sourced relative to this folder. 

- Currently the data needs to be supplied as a .rds file in the "data"-folder. The data needs to be in long format. The individual id variable should be numeric and start at 1 without gaps, i.e. 1:N. The alternative and choice task variables should also start indexing at 1, i.e 1:J and 1:T. See "data.demo.rds" for an example. If you want to use data stored in other formats you will need to change the corresponding line in the fn.setup.data() in the "r.script.methods.data.R"" -file. 

- You specify your model by changing the options in "r.script.run.eclc.mnl.R". This is the only file you need to change to run your model. The last line of this file: "fn.run.model()", will initiate the run sequence. 

- A completed model will store the model object as a .rds file and the output as a .txt file for easy inspection and use in post-estimation.

################################################################################
##  References
################################################################################

Hole, A. R., 2011, A discrete choice model with endogeneous attribute attendance, Economics Letters, 110(3):203-205


