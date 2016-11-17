setwd('/home/dcries')
Rcpp::sourceCpp('ebmodel/ppred_analysis.cpp')
Rcpp::sourceCpp('ebmodel/bivar_fullmcmc.cpp')
source('ebmodel/base_fcn.R')
source("ebmodel/dpmm_simulation.R")

dpmm_2_3 <- run_dpmmsim(2,3)

save(dpmm_2_3,file="dpmm_2_3.RData")
