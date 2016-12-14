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
#params <- c(100,50,300,14,-7,-200,8,-5)
simdata <- generate_data4(params,dist=1,nrep=2)
#simdata2 <- generate_data2(params,dist=1)

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
nreps <- 10000
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



ck <- 0.4
initial <- list(currentkee=currentkee,currentkes=currentkes,ck=ck,knotsee=knotsee,
                knotses=knotses,currentxee=currentxee,currentxes=currentxes,
                currentpredee=currentpredee,currentpredes=currentpredes,
                currentsigma2ee=currentsigma2ee,
                currentsigma2es=currentsigma2es,
                currentbetaee=currentbetaee,currentbetaes=currentbetaes)



prior <- list(lambda=lambda,ae=ae,be=be,a_alp=a_alp,
              b_alp=b_alp,d=d,m=m,v2=v2,psi=psi,Vb=Vb,Mb=Mb)


chain1=mcmc_naive(yee,yes,Z,initial,prior,nreps,burn,maxkt,my_bs,my_qp)
