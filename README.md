# GitHub
A spatio-temporal approach for radar image analysis
Data: March 14, 2017

#******************************************************
#******************************************************
The R code is based on the following paper submitted to the Annals of Applied Statistics:
"A Spatio-Temporal Modeling Approach for Weather Radar Reflectivity Data and Its Applications"
Author1: Xiao Liu (IBM Thomas J. Watson Research Center) 
Author2: Vikneswaran Gopal (Department of Statistics and Applied Probability, National University of Singapore)
Author3: Jayant Kalagnanam (IBM Thomas J. Watson Research Center) 
#******************************************************
#******************************************************

#******************************************************
#******************************************************
# Contact information: Xiao Liu, liuxiao314923@gmail.com
#******************************************************
#******************************************************

#******************************************************
#******************************************************
The main code:
“new_QPF_IBM_v1.R”

Sub-rountines:
“new_QPF_IBM_v1_Kernel.R”
“SubRoutine.R”
“TREC_function.R”
“TREC_function_P1.R”
“TREC_function_P2.R”
“new_QPF_IBM_v1_pred.R”
#******************************************************
#******************************************************

#******************************************************
#******************************************************
R packages needed:
sp, rgdal, rgeos, gstat, RColorBrewer, spacetime, fields, timeSeries, MASS, snow, mvtnorm, matrixcalc

A special package: selextR
A special packages created by Vik Gopal, which can be installed from "selextR_0.1.tar.gz"
#******************************************************
#******************************************************

#******************************************************
#******************************************************
Sample Data
A sample dataset is provided in the folder: data
#******************************************************
#******************************************************

#******************************************************
#******************************************************
Output
All output are saved in the folder: output
#******************************************************
#******************************************************

