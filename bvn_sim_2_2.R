setwd('/home/dcries')
Rcpp::sourceCpp('ebmodel/ppred_analysis.cpp')
Rcpp::sourceCpp('ebmodel/bivar_bvnmcmc1_qp.cpp')
source('ebmodel/base_fcn.R')
source("ebmodel/bvn_simulation.R")

bvn_2_2 <- run_bvnsim(2,2)

save(bvn_2_2,file="bvn_2_2.RData")
