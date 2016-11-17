setwd('/home/dcries')
Rcpp::sourceCpp('ebmodel/ppred_analysis.cpp')
Rcpp::sourceCpp('ebmodel/bivar_fullmcmc.cpp')
source('ebmodel/base_fcn.R')
source("ebmodel/dpmm_simulation.R")

dpmm_4_3 <- run_dpmmsim(4,3)

save(dpmm_4_3,file="dpmm_4_3.RData")
