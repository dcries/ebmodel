setwd('/home/dcries')
Rcpp::sourceCpp('ebmodel/ppred_analysis.cpp')
Rcpp::sourceCpp('ebmodel/bivar_bvnmcmc1_qp.cpp')
source('ebmodel/base_fcn.R')
source("ebmodel/bvn_simulation.R")

bvn_2_4 <- run_bvnsim(4,2)

save(bvn_2_4,file="bvn_4_2.RData")
