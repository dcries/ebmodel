run_bvnsim <- function(nreplicates,dist){
  library(MASS)
  #library(ggplot2)
  #library(gridExtra)
  library(mgcv) 
  library(reshape)
  library(mvtnorm)
  library(MCMCpack)
  library(splines)
  library(fGarch)
  #library(xtable)
  library(quadprog)
  
  #Rcpp::sourceCpp('C:/Users/dcries/github/ebmodel/ppred_analysis.cpp')
  #source('C:\\Users\\dcries\\github\\ebmodel\\base_fcn.R')
  # setwd('/home/dcries')
  # Rcpp::sourceCpp('/ebmodel/ppred_analysis.cpp')
  # Rcpp::sourceCpp('/ebmodel/bivar_bvnmcmc1_qp.cpp')
  # source('/ebmodel/base_fcn.R')
  #
  # !!!!!!!!!!!!!!!!!!
  # this value is hardcoded here
  n <- 300
  # !!!!!!!!!!!!!!!!!!
  
  set.seed(936)
  #number of simulated data sets
  nsim <- 200
  #number of posterior predictive iterations
  msim <- 1000
  #burn in
  nadapt <- 2000
  #iterations
  niter <- 10000
  #n chains
  nchain <- 1
  #pick which individuals interested in
  ind <- c(72,188,259)
  nind <- length(ind)
  
  pmse <- matrix(0,nrow=2,ncol=nsim)
  model <- matrix(0,nrow=15,ncol=nsim)
  
  indcheck <- matrix(0,nrow=nsim,ncol=4*length(ind))
  
  params <- c(100,50,300,14,-7,-200,8,-5)
  #names for X diagnostics
  wnames <- c("EE Min","EE 5%tile","EE 10%tile","EE 25%tile","EE 75%tile","EE 90%tile","EE 95%tile","EE Max",
              "ES Min","ES 5%tile","ES 10%tile","ES 25%tile","ES 75%tile","ES 90%tile","ES 95%tile","ES Max")
  
  #number of mcmc iterations after burn
  ureps <- niter
  #tuning burnin
  burn <- nadapt
  #number of iterations needed
  nreps <- ureps+burn 
  #inital number of knots
  currentkee <- 2
  currentkes <- 2
  #inital knot locations
  #knots <- sort(x)[c(50,125,200,250)]
  #specified by Denison
  ck <- 0.4
  #number of continuous derivatives -1
  l <- 3
  #number of components
  h <- 1#10
  #sd for random walk
  maxkt <- 15
  
  
  #current latent variables x, used in lm(y~bs(x))
  currentmuee <- rep(2600,h)
  currentmues <- rep(0,h)
  #currentmuee <- 2600#runif(h,1800,3800)
  #currentmues <- 0#runif(h,-400,400)
  currentpi <- rep(1/h,h)
  currentzeta <- sample(1:h,n,replace=T)
  currentv <- rep(0.3,h)

  currentsigma2ee <- 400^2
  currentsigma2ve <- 250^2
  currentsigma2es <- 220^2
  currentsigma2vs <- 80^2
  currentalpha <- 1
  tunevar = c(100^2,25^2)#c(100^2,25^2)
  tunecor = -0.2
  currentsigma2x <- c(400^2,400^2)
  currentbetaee <- matrix(c(0,0,0),nrow=1)
  currentbetaes <- matrix(c(0,0,0),nrow=1)
  
  
  #priors
  #lambda is value for mean of prior for number of knots
  lambda <- 1
  #priors for sigmae
  ae <- 0.01;be <- 0.01
  #priors for sigmav
  av <- 0.01; bv <- 0.01
  m <- c(2400,0)
  v2 <- matrix(c(3000^2,0,0,1000^2),ncol=2,byrow=FALSE)
  a_alp <- 1
  b_alp <- 1
  d <- 3#3
  psi <- diag(2)
  #prior variance for all coefficients
  Vb <- 100000
  Mb <- 0
  
  
  ck <- 0.4
  
  prior <- list(lambda=lambda,ae=ae,be=be,av=av,bv=bv,a_alp=a_alp,
                b_alp=b_alp,d=d,m=m,v2=v2,psi=psi,Vb=Vb,Mb=Mb)
  
  
  for(i in 1:nsim){
    params <- c(100,50,300,14,-7,-200,8,-5)
    #params <- c(100,50,300,14,-7,-200,8,-5)
    simdata <- generate_data5(params,dist=dist,nrep=nreplicates)
    
    yee <- simdata$yee
    yes <- simdata$yes
    yeeb <- rowMeans(yee)
    yesb <- rowMeans(yes)
    xee <- simdata$xee
    xes <- simdata$xes
    wee <- simdata$wee
    wes <- simdata$wes
    n <- length(yeeb)
    
    zg <- simdata$zg #<- rbinom(n,1,0.5) #gender indicator
    zb <- simdata$zb#<- rnorm(n,27,5) #bmi
    za <- simdata$za#<- runif(n,20,40) #age
    
    Z= cbind(zg,zb,za)
    
    knotsee <- sort(wee)[c(50,200)]
    knotses <- sort(wes)[c(50,200)]
    
    currentpredee <- wee[,1]
    currentpredes <- wes[,1]
    
    currentxee <- wee[,1]
    currentxes <- wes[,1]
    
    
    initial <- list(currentkee=currentkee,currentkes=currentkes,ck=ck,knotsee=knotsee,
                    knotses=knotses,currentxee=currentxee,currentxes=currentxes,
                    currentv=currentv,currentpi=currentpi,currentalpha=currentalpha,
                    currentzeta=currentzeta,currentpredee=currentpredee,currentpredes=currentpredes,
                    currentmuee=currentmuee,currentmues=currentmues,currentsigma2ee=currentsigma2ee,
                    currentsigma2es=currentsigma2es,currentsigma2ve=currentsigma2ve,
                    currentsigma2vs=currentsigma2vs,currentsigma2x=currentsigma2x,
                    tunevar=tunevar,currentbetaee=currentbetaee,currentbetaes=currentbetaes,
                    tunecor=tunecor)
    
    sample=mcmc_bvn_qp(yee,yes,wee,wes,Z,initial,prior,nreps,burn,maxkt,my_bs,my_qp)
    
    pp_indicies <- sample(1:niter,msim)
    
    sample_r <- list(meanfcnee=sample$meanfcnee[pp_indicies,],meanfcnes=sample$meanfcnes[pp_indicies,],
                     betaee=sample$betaee[pp_indicies,],betaes=sample$betaes[pp_indicies,])
    
    pp <- pp_pmse(sample_r,yeeb,yesb,Z)
    
    pmse[,i] <- colMeans(pp)
    model[,i] <- c(mean(sample$muee),mean(sample$mues),sqrt(mean(sample$sigma2eee)),
                   sqrt(mean(sample$sigma2ees)),sqrt(mean(sample$sigma2vee)),
                   sqrt(mean(sample$sigma2ves)),sqrt(mean(sample$sigma2xee)),
                   sqrt(mean(sample$sigma2xes)),mean(sample$corrx),
                   colMeans(sample$betaee),colMeans(sample$betaes))
    
    indcheck[i,1] <- simdata$xee[ind[1]]
    indcheck[i,2] <- mean(sample$latentxee[,ind[1]])
    indcheck[i,3] <- simdata$xes[ind[1]]
    indcheck[i,4] <- mean(sample$latentxes[,ind[1]])
    
    indcheck[i,5] <- simdata$xee[ind[2]]
    indcheck[i,6] <- mean(sample$latentxee[,ind[2]])
    indcheck[i,7] <- simdata$xes[ind[2]]
    indcheck[i,8] <- mean(sample$latentxes[,ind[2]])
    
    indcheck[i,9] <- simdata$xee[ind[3]]
    indcheck[i,10] <- mean(sample$latentxee[,ind[3]])
    indcheck[i,11] <- simdata$xes[ind[3]]
    indcheck[i,12] <- mean(sample$latentxes[,ind[3]])
    
    print(i)
    
  }
  
  pmse <- t(pmse)
  model <- t(model)
  
  return(list(pmse=pmse,model=model,indcheck=indcheck))
  
}


