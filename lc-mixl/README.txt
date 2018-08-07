################################################################################
##  Name:         README.txt
##  Created:      2018.08.02
##  Last edited:  2018.08.07
##  Author:       Erlend Dancke Sandorf
################################################################################

################################################################################
Generic:

- Keep folder structure the same

- You should only need to change the file: r.script.run.lc.mixl.R
  This is where you change the options of estimation and specify the data to use

- The data needs to be saved as .rds and in long format. This can also be changed in the fn.setup.data()

- The parameter order is as follows: fixed params, means class1, means class 2, ..., sds class1, sds class2, ..., 

- The model allows the parameter distributions to be correlated (within-class) or uncorrelated and interactions with the means to be either generic or class specific

################################################################################
References:

Greene, W. & Hensher, D. A., 2013, Revealing additional dimensions of preference heterogeneity in a latent class mixed multinomial logit model, Applied Economics, 45(14):1897-1902
