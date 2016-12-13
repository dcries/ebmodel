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
Rcpp::sourceCpp('C:\\Users\\dcries\\github\\ebmodel\\ppred_analysis.cpp')
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

#number of mcmc iterations after burn
ureps <- 2000
#tuning burnin
burn <- 00
#number of iterations needed
nreps <- ureps+burn 
#inital number of knots
currentkee <- 2;currentkee2 <- 1;currentkee3 <- 5
currentkes <- 2;currentkes2 <- 5;currentkes3 <- 1
#inital knot locations
#knots <- sort(x)[c(50,125,200,250)]
knotsee <- sort(wee)[c(50,200)];knotsee2 <- sort(wee)[c(200)];knotsee3 <- sort(wee)[c(10,50,100,200,250)]
knotses <- sort(wes)[c(50,200)];knotses2 <- sort(wes)[c(10,50,100,200,250)];knotses3 <- sort(wes)[c(50)]
#specified by Denison
ck <- 0.4
#number of continuous derivatives -1
l <- 3
#number of components
h <- 1#10
#sd for random walk
maxkt <- 15


#current latent variables x, used in lm(y~bs(x))
currentxee <- rowMeans(wee);currentxee2 <- yee[,2];currentxee3 <- wee[,2]+rnorm(n,0,300)
currentxes <- rowMeans(wes);currentxes2 <- yes[,1];currentxes3 <- wes[,2]+rnorm(n,0,100)
#currentmuee <- 2600#runif(h,1800,3800)
#currentmues <- 0#runif(h,-400,400)
currentpredee <- wee[,1]
currentpredes <- wes[,1]
currentsigma2ee <- 400^2;currentsigma2ee2 <- 200^2;currentsigma2ee3 <- 600^2
currentsigma2es <- 220^2;currentsigma2es2 <- 120^2;currentsigma2es3 <- 320^2
currentbetaee <- matrix(c(0,0,0),nrow=1);currentbetaee2 <- matrix(c(-50,-15,-15),nrow=1);currentbetaee3 <- matrix(c(100,15,15),nrow=1)
currentbetaes <- matrix(c(0,0,0),nrow=1);currentbetaes2 <- matrix(c(-100,-15,-15),nrow=1);currentbetaes3 <- matrix(c(100,15,15),nrow=1)


#priors
#lambda is value for mean of prior for number of knots
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
