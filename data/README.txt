################################################################################
##  Name:         README.txt
##  Created:      2018.08.31
##  Last edited:  2018.08.31
##  Author:       Erlend Dancke Sandorf
################################################################################

################################################################################
##  Data structure
################################################################################
The data needs to be stored in long format (one row per alternative). The data must include an id-variable for the individual that is continuous starting at 1, a choice task variable that is continuous starting at 1 and an alternative variable that is continuous starting at 1. 

################################################################################
##  Cold water coral data
################################################################################
Dataset containing choices from 397 individuals answering 12 choice tasks containing 2 experimentally designed alternatives and a status quo. The data is "incomplete", i.e. not all individuals answered all choice tasks. There are some missing values for some of the socio-demographic variables. 

Individuals: 	397
Choice tasks: 	12
Alternatives: 	3
Observations:	4683
Variables:	26

id		- Individual identification variable
ct		- Choice task identification variable
alt		- Alternative identification variable
choice		- Choice indicator
size		- Size of the protected area (Attribute 1)
small		- Indicator for a small increase (Attribute 1)
large		- Indicator for a large increase (Attribute 1)
oil		- Indicator for whether the proposed protected area is important for the oil and gas industry (Attribute 2)
fish		- Indicator for whether the proposed protected area is important for the fishing industry (Attribute 2)
hab		- Indicator for whether the proposed area is important habitat for fish (Attribute 3)
cost		- Cost measured as a lump sum increase in annual taxes (Attribute 4)
complete	- Indicates the number of completed choice tasks
out		- Indicator for the choice tasks where no choice was made
treatment	- Indicator for whether a respondent received his or her score
score		- Sce on the CWC quiz
male		- Indicator for male respondent

################################################################################
##  References
################################################################################
Data:
Aanesen, M., 2018, Cold-water coral protection in Norway: A choice experiment, DataverseNO, Version 1, https://doi.org/10.18710/Q39V7O

Original publication: 
Aanesen, M., Armstrong, C., Czajkowski, M., Falk-Petersen, J., Hanley, N. & Navrud, S., 2015, Willingness to pay for unfamiliar public goods: Preserving cold-water coral in Norway, Ecological Economics, 112:53-67






