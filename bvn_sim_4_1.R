setwd('/home/dcries')
Rcpp::sourceCpp('ebmodel/ppred_analysis.cpp')
Rcpp::sourceCpp('ebmodel/bivar_bvnmcmc1_qp.cpp')
source('ebmodel/base_fcn.R')
source("ebmodel/bvn_simulation.R")

bvn_4_1 <- run_bvnsim(4,1)

save(bvn_4_1,file="bvn_4_1.RData")
