library(MASS)
library(ggplot2)
library(gridExtra)
library(mgcv) 
library(reshape)
library(mvtnorm)
library(MCMCpack)
library(splines)
library(fGarch)
library(xtable)
library(quadprog)

#source('//my.files.iastate.edu/Users/dcries/Desktop/research/rprograms/base_fcn.R', echo=FALSE)#Rcpp::sourceCpp('//my.files.iastate.edu/Users/dcries/Desktop/research/rprograms/trial_and_simulated/bivar_fullmcmc2.cpp')
source('C:\\Users\\dcries\\github\\ebmodel\\base_fcn.R')
Rcpp::sourceCpp('C:\\Users\\dcries\\github\\ebmodel\\bivar_naivemcmc.cpp')

params <- c(100,50,300,14,-7,-200,8,-5)
simdata <- generate_data5(params,dist=1,nrep=2)

zg <- simdata$zg #<- rbinom(n,1,0.5) #gender indicator
zb <- simdata$zb#<- rnorm(n,27,5) #bmi
za <- simdata$za#<- runif(n,20,40) #age

Z= cbind(zg,zb,za)


yeeb <- rowMeans(yee)
xee <- rowMeans(wee)
yesb <- rowMeans(yes)
xes <- rowMeans(wes)
np <- ncol(Z)
nr <- ncol(yee)
lambda <- 1
#priors for sigmae
ae <- 0.01;be <- 0.01

a_alp <- 1
b_alp <- 1
d <- 3#3
psi <- diag(2)
#prior variance for all coefficients
Vb <- 100000
Mb <- 0


pick_indicies <- function(l0,full,sel){
  l1 <- l0+1
  temp <- match(sel,full) + matrix(rep(-l1:l1,each=length(sel)),ncol=2*l1+1,byrow=FALSE)
  out <- c(as.numeric(temp),1:l0,(length(full)-l0):length(full))
  return(out[out>0])
}

log_q <- function(pmean,sigma2,y){
  n <- length(y)
  return(-n*log(sqrt(sigma2)) - 1/(2*sigma2)*sum((y-pmean)^2))
}

#number of mcmc iterations
nreps <- 2000
#inital number of knots
kes <- rep(0,nreps)
kee <- rep(0,nreps)

currentkee <- 2
currentkes <- 2

sigma2ee <- rep(0,nreps)
sigma2es <- rep(0,nreps)

currentsigma2ee <- 400
currentsigma2es <- 200

#inital knot locations
knotsee <- sort(xee)[c(100,200)]
knotses <- sort(xes)[c(100,200)]

#lambda is value for mean of prior for number of knots
lambda <- 2
#priors for sigma
a <- 0.01;b <- 0.01
#specified by Denison
ck <- 0.4
#number of continuous derivatives -1
l <- 3
#current predictions for fucntinos for likelihood ratio
currentpredee <- rowMeans(wee)
currentbetaee <- rep(0,ncol(Z))
#mean fcn
meanfcnee <- matrix(0,nrow=n,ncol=nreps)
#prediction locations to check mixing
betaee <- matrix(0,nrow=nreps,ncol=ncol(Z))
currentpredes <- rowMeans(wes)
currentbetaes <- rep(0,ncol(Z))
#mean fcn
meanfcnes <- matrix(0,nrow=n,ncol=nreps)
#prediction locations to check mixing
betaes <- matrix(0,nrow=nreps,ncol=ncol(Z))

for(i in 1:nreps){
  
  bk <- ck*min(c(1,dpois(currentkee+1,lambda)/dpois(currentkee,lambda)))
  dk <- ck*min(c(1,dpois(currentkee,lambda)/dpois(currentkee+1,lambda)))
  unavailable <- pick_indicies(l,sort(xee),sort(knotsee))
  u <- runif(1)
  if(u <= bk){ #birth step
    new <- sample(sort(xee)[-unavailable],1)
    knotspropee <- sort(c(knotsee,new))
    #bmatrix <- bs(xee,knots=knotspropee)
    bmatrix <- my_bs(xee,knots=knotspropee)
    
    #m1 <- lm(yeeb~bmatrix)
    m1 <- my_lm_qp(yeeb,bmatrix,my_qp)
    
    acceptprob <- log_q(unlist(m1["preds"]),currentsigma2ee,yeeb) - log_q(currentpredee,currentsigma2ee,yeeb) +
      log(n - (2*(l+1) + (currentkee)*(2*l+1))) - log(n)#- 0.5*log(n) 
    #acceptprob <- log_q(predict(m1),currentsigma2ee,yeeb) - log_q(currentpredee,currentsigma2ee,yeeb) +
    #  log(n - (2*(l+1) + (currentkee)*(2*l+1))) - log(n)#- 0.5*log(n) 
    
    if(log(runif(1)) < acceptprob){
      currentpredee <- unlist(m1["preds"])
      #currentpredee <- predict(m1)
      knotsee <- knotspropee
    }
    
  } 
  else if(u <= bk+dk & length(knotsee)>1){ #death step
    rmv <- sample(knotsee,1)
    knotspropee <- sort(knotsee[-which(knotsee==rmv)])
    #bmatrix <- bs(xee,knots=knotspropee)
    bmatrix <- my_bs(xee,knots=knotspropee)
    
    #m1 <- lm(yeeb~bmatrix)
    m1 <- my_lm_qp(yeeb,bmatrix,my_qp)
    
     acceptprob <- log_q(unlist(m1["preds"]),currentsigma2ee,yeeb) - log_q(currentpredee,currentsigma2ee,yeeb) -
       log(n - (2*(l+1) + currentkee*(2*l+1))) + log(n)#- 0.5*log(n)
    #acceptprob <- log_q(predict(m1),currentsigma2ee,yeeb) - log_q(currentpredee,currentsigma2ee,yeeb) -
    #  log(n - (2*(l+1) + currentkee*(2*l+1))) + log(n)#- 0.5*log(n)
    
    if(log(runif(1)) < acceptprob){
      currentpredee <- unlist(m1["preds"])
      #currentpredee <- predict(m1)
      
      knotsee <- knotspropee
    }
    
  }
  else{ #move step
    chg <- sample(knotsee,1)
    knotspropee <- knotsee
    knotspropee[which(knotspropee==chg)] <- sample(sort(xee)[-unavailable],1)
    knotspropee <- sort(knotspropee)
    #bmatrix <- bs(xee,knots=knotspropee)
    bmatrix <- my_bs(xee,knots=knotspropee)
    
    #m1 <- lm(yeeb~bmatrix)
    m1 <- my_lm_qp(yeeb,bmatrix,my_qp)
    
    acceptprob <- log_q(unlist(m1["preds"]),currentsigma2ee,yeeb) - log_q(currentpredee,currentsigma2ee,yeeb)#- 0.5*log(n)
    #acceptprob <- log_q(predict(m1),currentsigma2ee,yeeb) - log_q(currentpredee,currentsigma2ee,yeeb)#- 0.5*log(n)
    
    if(log(runif(1)) < acceptprob){
      currentpredee <- unlist(m1["preds"])
      #currentpredee <- predict(m1)
      
      knotsee <- knotspropee
    }
  }
  
  
  #es spline step
  bk <- ck*min(c(1,dpois(currentkes+1,lambda)/dpois(currentkes,lambda)))
  dk <- ck*min(c(1,dpois(currentkes,lambda)/dpois(currentkes+1,lambda)))
  unavailable <- pick_indicies(l,sort(xes),sort(knotses))
  u <- runif(1)
  if(u <= bk){ #birth step
    new <- sample(sort(xes)[-unavailable],1)
    knotspropes <- sort(c(knotses,new))
    #bmatrix <- bs(xes,knots=knotspropes)
    bmatrix <- my_bs(xes,knots=knotspropes)
    
    #m1 <- lm(yesb~bmatrix)
    m1 <- my_lm_qp(yesb,bmatrix,my_qp)
    
     acceptprob <- log_q(unlist(m1["preds"]),currentsigma2es,yesb) - log_q(currentpredes,currentsigma2es,yesb) +
       log(n - (2*(l+1) + (currentkes)*(2*l+1))) - log(n)#- 0.5*log(n)
    #acceptprob <- log_q(predict(m1),currentsigma2es,yesb) - log_q(currentpredes,currentsigma2es,yesb) +
    #  log(n - (2*(l+1) + (currentkes)*(2*l+1))) - log(n)#- 0.5*log(n)
    
   if(log(runif(1)) < acceptprob){
      currentpredes <- unlist(m1["preds"])
      #currentpredes <- predict(m1)
      
      knotses <- knotspropes
    }
    
  }
  else if(u <= bk+dk & length(knotses)>1){ #death step
    rmv <- sample(knotses,1)
    knotspropes <- sort(knotses[-which(knotses==rmv)])
    #bmatrix <- bs(xes,knots=knotspropes)
    bmatrix <- my_bs(xes,knots=knotspropes)
    
    #m1 <- lm(yesb~bmatrix)
    m1 <- my_lm_qp(yesb,bmatrix,my_qp)
    
     acceptprob <- log_q(unlist(m1["preds"]),currentsigma2es,yesb) - log_q(currentpredes,currentsigma2es,yesb) -
       log(n - (2*(l+1) + currentkes*(2*l+1))) + log(n)#- 0.5*log(n)
    #acceptprob <- log_q(predict(m1),currentsigma2es,yesb) - log_q(currentpredes,currentsigma2es,yesb) -
    #  log(n - (2*(l+1) + currentkes*(2*l+1))) + log(n)#- 0.5*log(n)
    
    if(log(runif(1)) < acceptprob){
      currentpredes <- unlist(m1["preds"])
      #currentpredes <- predict(m1)
      
      knotses <- knotspropes
    }
    
  }
  else{ #move step
    chg <- sample(knotses,1)
    knotspropes <- knotses
    knotspropes[which(knotspropes==chg)] <- sample(sort(xes)[-unavailable],1)
    knotspropes <- sort(knotspropes)
    #bmatrix <- bs(xes,knots=knotspropes)
    bmatrix <- my_bs(xes,knots=knotspropes)
    
    #m1 <- lm(yesb~bmatrix)
    m1 <- my_lm_qp(yesb,bmatrix,my_qp)
    
    acceptprob <- log_q(unlist(m1["preds"]),currentsigma2es,yesb) - log_q(currentpredes,currentsigma2es,yesb)#- 0.5*log(n)
    #acceptprob <- log_q(predict(m1),currentsigma2es,yesb) - log_q(currentpredes,currentsigma2es,yesb)#- 0.5*log(n)
    
    if(log(runif(1)) < acceptprob){
      currentpredes <- unlist(m1["preds"])
      #currentpredes <- predict(m1)
      
      knotses <- knotspropes
    }
  }
  
  
  
  yeebdiff = yeeb - currentpredee
  yesbdiff = yesb - currentpredes
  
  Vee = solve(t(Z)%*%Z + (1/Vb)*diag(np))
  Mee = Vee%*%((1/Vb)*Mb + t(Z)%*%yeebdiff)
  Ves = solve(t(Z)%*%Z + (1/Vb)*diag(np))
  Mes = Ves%*%((1/Vb)*Mb + t(Z)%*%yesbdiff)
  e2ee = 0
  e2es = 0
  
  for(ii in 1:nr){
    #e2ee = e2ee + sum((yee[,ii] - currentpredee - Z%*%currentbetaee)^2)/2.0
    e2ee = e2ee + t(yee[,ii]- currentpredee)%*%(yee[,ii] - currentpredee)
    e2es = e2es + sum((yes[,ii] - currentpredes - Z%*%currentbetaes)^2)/2.0
  }

  currentsigma2ee = rinvgamma(1,ae+nr*n/2.0,(be+(e2ee- t(Mee)%*%solve(Vee)%*%Mee)/2))
  currentsigma2es = rinvgamma(1,ae+nr*n/2.0,(be+e2es))
  

  currentbetaee = mvrnorm(1,Mee,currentsigma2ee*Vee)
  currentbetaes = mvrnorm(1,Mes,currentsigma2es*Ves)
  
  #currentsigma2 <- rinvgamma(1,a+n/2,b+sum((y-currentpred)^2)/2)
  
  kee[i] <- length(knotsee) #number of knots in model 
  meanfcnee[,i] <- currentpredee 
  sigma2ee[i] <- currentsigma2ee
  betaee[i,] <- currentbetaee
  currentkee <- kee[i]
  
  kes[i] <- length(knotses) #number of knots in model
  meanfcnes[,i] <- currentpredes
  sigma2es[i] <- currentsigma2es
  betaes[i,] <- currentbetaes
  currentkes <- kes[i]
  
  #print(knots)
  
}
