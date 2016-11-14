#Fit multiple different models from realdata_sim and calc DIC, posterior pred
#measures to assess 

#save.image("//my.files.iastate.edu/Users/dcries/Desktop/jags_mcmc1.RData")

library(rstan)
library(rjags)
library(MASS)
library(ggplot2)
library(gridExtra)
library(mgcv) 
library(reshape)
library(mvtnorm)
library(fGarch)

setwd("U:/Desktop/research/data")
#source("U:\\Desktop\\research\\rprograms\\base_fcn.R")
#Rcpp::sourceCpp('//my.files.iastate.edu/Users/dcries/Desktop/research/rprograms/trial_and_simulated/ppred_analysis.cpp')
#source('\\\\my.files.iastate.edu\\Users\\dcries\\Desktop\\research\\rprograms\\base_fcn.R')
Rcpp::sourceCpp('C:\\Users\\dcries\\github\\ebmodel\\ppred_analysis.cpp')
source('C:\\Users\\dcries\\github\\ebmodel\\base_fcn.R')
eb <- read.csv("ebs_2_23.csv")

params <- c(100,50,300,14,-7,-200,8,-5)
simdata <- generate_data4(params,dist=1)
simdata$d <- 3 # for inv-wish prior
simdata$I <- diag(2) #for inv-wish prior
simdata$n <- length(simdata$xee)
simdata$H <- 20 #number of clusters for DP
simdata$alpha <- 2
simdata$knots <- seq(1400,4000,length=6)
simdata$knots2 <- seq(-800,800,length=6)
simdata$K <- length(simdata$knots) #number of knots


mu0 <- matrix(c(2800,0,2400,0,2000,0),ncol=2,byrow=TRUE)
mucov <- array(0,dim=c(2,2,3))
diag(mucov[,,1]) <- c(1/(200^2),1/(200^2))
diag(mucov[,,2]) <- c(1/(200^2),1/(200^2))
diag(mucov[,,3]) <- c(1/(200^2),1/(200^2))

mu0dp <- c(2400,0)
mucovdp <- matrix(c(1/(200^2),0,0,1/(200^2)),ncol=2,byrow=TRUE)

simdata$mu0 <- mu0
simdata$mucov <- mucov
simdata$a <- rep(1,2)
simdata$mu0dp <- mu0dp
simdata$mucovdp <- mucovdp

#quantiles of xei, range of xee, P(xee < 2000), P(xes < 0), ie burning more than gaining,
#t_sum <- c(quantile(simdata$xei,probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)),range(simdata$xee),sum(simdata$xee < 2000)/length(simdata$xee),sum(simdata$xes < 0)/length(simdata$xes),range(simdata$xee+simdata$xei))
t_sum <- c(min(simdata$xei),quantile(simdata$xei,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$xei),min(simdata$xee),quantile(simdata$xee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$xee))
t_sum2 <- c(sum(simdata$xee < 2000)/length(simdata$xee),sum(simdata$xes < 0)/length(simdata$xes),range(simdata$xee+simdata$xei),quantile(simdata$xee,probs=c(0.75))-quantile(simdata$xee,probs=c(0.25)),quantile(simdata$xei,probs=c(0.75))-quantile(simdata$xei,probs=c(0.25)),range(simdata$xes))
#t_sumee <- c(quantile(simdata$xee,probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)),range(simdata$xei),quantile(simdata$xee,probs=c(0.75))-quantile(simdata$xee,probs=c(0.25)),quantile(simdata$xei,probs=c(0.75))-quantile(simdata$xei,probs=c(0.25)),range(simdata$xes))
#discrepancy measures for W measurements
wsum <- c(min(simdata$wee),quantile(simdata$wee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$wee),min(simdata$wes),quantile(simdata$wes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$wes))
#discrepancy meaures for Y measurements
ysum <- c(min(simdata$yee),quantile(simdata$yee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$yee),min(simdata$yes),quantile(simdata$yes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$yes))


#names for X diagnostics
xnames <- c("EI Min","EI 5%tile","EI 10%tile","EI 25%tile","EI 75%tile","EI 90%tile","EI 95%tile","EI Max",
            "EE Min","EE 5%tile","EE 10%tile","EE 25%tile","EE 75%tile","EE 90%tile","EE 95%tile","EE Max")

wnames <- c("EE Min","EE 5%tile","EE 10%tile","EE 25%tile","EE 75%tile","EE 90%tile","EE 95%tile","EE Max",
            "ES Min","ES 5%tile","ES 10%tile","ES 25%tile","ES 75%tile","ES 90%tile","ES 95%tile","ES Max")

#covariate values for pmse for checking regression model fits
set.seed(114)
#known covariates
n <- length(simdata$xee)
zg <- simdata$zg #<- rbinom(n,1,0.5) #gender indicator
zb <- simdata$zb#<- rnorm(n,27,5) #bmi
za <- simdata$za#<- runif(n,20,40) #age
xee <- simdata$xee
xes <- simdata$xes

weeb <- rowMeans(simdata$wee)
wesb <- rowMeans(simdata$wes)
yeeb <- rowMeans(simdata$yee)
yesb <- rowMeans(simdata$yes)

#need "truth" for sigmas for Y and W
#for dist=1, ie Normal
#sigmauee, .1*xee + delta ee?, with both sd(wee1-xee) = 286, sd(wee2-xee) = 271
#sigamues .1*xes+100 + delta es? sd(wes1-xes) = 196, sd(wes2-xes)=205
#sigmaeee .14*xee + delta ee? sd(yee1-(be0+eecurve(xee)+geg*zg+geb*zb+gea*za))=401,388 
#sigmaees .25*xes+250 + delta es, sd(yes2-(bs0+escurve(xes)+gig*zg+gib*zb+gia*za))=320,349

#for dist=2, skewed normal
#sigmauee sd(wee1-xee) = 262,260
#sigmaues sd(wes1-xes) = 178,172
#sigmaeee sd() =382,376
#sigmaees sd = 303,310

#for dist=3, bimodal centered at 0
#sigmauee sd = 196
#sigmaues sd = 173
#sigmaeee sd = 312
#sigmaees sd = 326

#-----------------------------------------------#
#-----------------Model 0-----------------------#
#naive model, assume w are measurements of X without error

modelj0 <- "
model{
for(i in 1:n){
for(j in 1:2){
yee[i,j] ~ dnorm(be0 + be1*x[i,1] + geg*zg[i] + geb*zb[i] + gea*za[i],taueee)
yes[i,j] ~ dnorm(bs0 + bs1*x[i,2] + gig*zg[i] + gib*zb[i] + gia*za[i],tauees)
#wee[i,j] ~ dnorm(x[i,1], tauuee)
#wes[i,j] ~ dnorm(x[i,2], tauues)
}
#x[i,1:2] ~ dmnorm(mux[],taux[,])
}

#mux[1] ~ dnorm(2200,1/300000)
#mux[2] ~ dnorm(0,1/100000)
taueee ~ dgamma(0.1,0.1) #~ dt(0,1,1)T(0,) half cauchy
tauees ~ dgamma(0.1,0.1)
#tauuee <- (1/sigmauee)^2
#tauues <- (1/sigmaues)^2
#taux[1:2,1:2] ~ dwish(I,d) #d should be dimension + 1
#sigmaeee ~ dt(0,1,1)T(0,) #standard deviations
#sigmaees ~ dt(0,1,1)T(0,)
sigmaeee <- 1/sqrt(taueee)
sigmaees <- 1/sqrt(tauees)
#sigmauee ~ dt(0,1,1)T(0,)
#sigmaues ~ dt(0,1,1)T(0,)

#sigmax <- inverse(taux[,])
#sigmaxee <- sqrt(sigmax[1,1])
#sigmaxes <- sqrt(sigmax[2,2])
#corrx <- sigmax[1,2]/(sigmaxee*sigmaxes)

be0 ~ dnorm(0,1/100000) #dnorm(100,1/10000)
be1 ~ dnorm(1,1/10000) #dnorm(1,1/10000)
bs0 ~ dnorm(0,1/100000) #dnorm(0,1/10000)
bs1 ~ dnorm(1,1/10000) #dnorm(1,1/10000)

geg ~ dnorm(0,1/100000)
geb ~ dnorm(0,1/100000)
gea ~ dnorm(0,1/100000)
gig ~ dnorm(0, 1/100000)
gib ~ dnorm(0,1/100000)
gia ~ dnorm(0,1/100000)
}
"

naivedata <- list(yee=simdata$yee[,1:2],yes=simdata$yes[,1:2],x=matrix(c(weeb,wesb),ncol=2,byrow=FALSE),n=length(simdata$xee),I=diag(2),d=3,zg=zg,za=za,zb=zb)
jm0 <- jags.model(textConnection(modelj0),data=naivedata,n.adapt=2000,n.chains=3)
js0 <- coda.samples(jm0,c("mux","sigmaeee","sigmaees","be0","be1","bs0","bs1","geg","geb","gea","gig","gib","gia"),n.iter=20000)
#dic0 <- dic.samples(jm0,n.iter=2000)
#gelman.diag(js0[,c( "sigmaeee","sigmaees","be0","be1","bs0","bs1","geg","geb","gea","gig","gib","gia")],multivariate=FALSE)
#summary(js0[,c("sigmaeee","sigmaees","be0","be1","bs0","bs1","geg","geb","gea","gig","gib","gia")])

jdf0 <- as.data.frame(as.matrix(js0))
#names(jdf0)[2:3] <- c("mux0","mux2")

#check0 <- matrix(0,nrow=nrow(jdf0),ncol=length(t_sum))
#check0b <- matrix(0,nrow=nrow(jdf0),ncol=length(t_sum2))
#check0w <- matrix(0,nrow=nrow(jdf0),ncol=length(wsum))
check0y <- matrix(0,nrow=nrow(jdf0),ncol=length(wsum))
dist0 <- matrix(0,nrow=nrow(jdf0),ncol=2)
data0 <- matrix(c(simdata$xee,simdata$xes),ncol=2,byrow=FALSE)
#data <- matrix(0,ncol=nrow(jdf0))
pmse0 <- matrix(0,ncol=2,nrow=nrow(jdf0))
for(i in 1:nrow(jdf0)){
  #covmat <- matrix(c(jdf0$sigmaxee[i]^2,jdf0$sigmaxee[i]*jdf0$sigmaxes[i]*jdf0$corrx[i],jdf0$sigmaxee[i]*jdf0$sigmaxes[i]*jdf0$corrx[i],jdf0$sigmaxes[i]^2),ncol=2,byrow=TRUE)
  
  #data <- mvrnorm(length(simdata$xee),c(jdf0$mux0[i],jdf0$mux2[i]),covmat)
  datayee <- rnorm(length(simdata$xee),jdf0$be0[i]+jdf0$be1[i]*weeb,jdf0$sigmaeee[i])
  datayes <- rnorm(length(simdata$xee),jdf0$bs0[i]+jdf0$bs1[i]*wesb,jdf0$sigmaees[i])
  
  # check0[i,1] <- min(data[,1]+data[,2])
  # check0[i,2:7] <- quantile(data[,1]+data[,2],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
  # check0[i,8] <- max(data[,1]+data[,2])
  # check0[i,9] <- min(data[,1])
  # check0[i,10:15] <- quantile(data[,1],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
  # check0[i,16] <- max(data[,1])
  # 
  # check0b[i,1] <- sum(data[,1] < 2000)/nrow(data)
  # check0b[i,2] <- sum(data[,2] < 0)/nrow(data)
  # check0b[i,3:4] <- range(2*data[,1] + data[,2])
  # check0b[i,5] <- quantile(data[,1],probs=c(0.75)) - quantile(data[,1],probs=c(0.25))
  # check0b[i,6] <- quantile(data[,1]+data[,2],probs=c(0.75)) - quantile(data[,1]+data[,2],probs=c(0.25))
  # check0b[i,7:8] <- range(data[,2])
  
  check0y[i,1] <- min(datayee)
  check0y[i,2:7] <- quantile(datayee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
  check0y[i,8] <- max(datayee)
  check0y[i,9] <- min(datayes)
  check0y[i,10:15] <- quantile(datayes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
  check0y[i,16] <- max(datayes)
  
  #dist0[i,1] <- sum(sqrt((data[,1]-simdata$xee)^2+(data[,2]-simdata$xes)^2))
  #dist0[i,2] <- sum(sqrt((data[,1]-data0[,1])^2+(data[,2]-data0[,2])^2))
  #data0 <- data
  
  pmse0[i,1] <- mean((yeeb - jdf0$be0[i]-jdf0$be1[i]*weeb-jdf0$geg[i]*zg-jdf0$geb[i]*zb-jdf0$gea[i]*za)^2)
  pmse0[i,2] <- mean((yesb - jdf0$bs0[i]-jdf0$bs1[i]*wesb-jdf0$gig[i]*zg-jdf0$gib[i]*zb-jdf0$gia[i]*za)^2)
}

#t_sum0 <- c(quantile(simdata$xei,probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)),range(simdata$xee),sum(simdata$xee < 2000)/length(simdata$xee),sum(simdata$xes < 0)/length(simdata$xes),range(simdata$xee+simdata$xei))
names0 <- c("$\\beta_{0,ee}$","$\\beta_{1,ee}$","$\\beta_{0,es}$","$\\beta_{1,es}$",
           "$\\sigma_{yee}$","$\\sigma_{yes}$",
           "$\\gamma_{1,ee}$","$\\gamma_{2,ee}$","$\\gamma_{3,ee}$",
           "$\\gamma_{1,es}$","$\\gamma_{2,es}$","$\\gamma_{3,es}$")
sumtable0 <- data.frame(summary(js0[,c("be0","be1","bs0","bs1", "sigmaeee","sigmaees","geg","geb","gea","gig","gib","gia")])$quantile)
rownames(sumtable0) <- names0
print(xtable(sumtable0), sanitize.text.function=function(x){x})


#check0 <- data.frame(check0)
#ppred16(check0,t_sum,xnames) #min and 5th percentile
#ppred8(check0b,t_sum2)
#ppred16(check0y,ysum,wnames)



#-----------------------------------------------#
#---------- Model 1-----------------------------#
#basic model
#assumes normality throughout, constant variance, y is unbiased

# modelj1 <- "
# model{
# for(i in 1:n){
# for(j in 1:2){
# yee[i,j] ~ dnorm(x[i,1],taueee)
# yes[i,j] ~ dnorm(x[i,2],tauees)
# wee[i,j] ~ dnorm(x[i,1], tauuee)
# wes[i,j] ~ dnorm(x[i,2], tauues)
# }
# x[i,1:2] ~ dmnorm(mux[],taux[,])
# }
# 
# mux[1] ~ dnorm(2200,1/300000)
# mux[2] ~ dnorm(0,1/100000)
# taueee ~ dgamma(0.1,0.1)#<- (1/sigmaeee)^2 #~ dt(0,1,1)T(0,) half cauchy
# tauees ~ dgamma(0.1,0.1)#<- (1/sigmaees)^2
# tauuee ~ dgamma(0.1,0.1)#<- (1/sigmauee)^2
# tauues ~ dgamma(0.1,0.1)#<- (1/sigmaues)^2
# taux[1:2,1:2] ~ dwish(I,d) #d should be dimension + 1
# sigmaeee <- sqrt(taueee)#~ dt(0,1,1)T(0,) #standard deviations
# sigmaees <- sqrt(tauees)#~ dt(0,1,1)T(0,)
# sigmauee <- sqrt(tauuee)#~ dt(0,1,1)T(0,)
# sigmaues <- sqrt(tauues)#~ dt(0,1,1)T(0,)
# 
# sigmax <- inverse(taux[,])
# sigmaxee <- sqrt(sigmax[1,1])
# sigmaxes <- sqrt(sigmax[2,2])
# corrx <- sigmax[1,2]/(sigmaxee*sigmaxes)
# }
# "
# 
# jm1 <- jags.model(textConnection(modelj1),data=simdata,n.adapt=2000,n.chains=3)
# js1 <- coda.samples(jm1,c("mux","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","x"),n.iter=20000)
# dic1 <- dic.samples(jm1,n.iter=2000)
# gelman.diag(js1[,c("mux[1]","mux[2]", "sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx")],multivariate=FALSE)
# summary(js1[,c("mux[1]","mux[2]","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx")])
# 
# jdf1 <- as.data.frame(as.matrix(js1))
# names(jdf1)[2:3] <- c("mux1","mux2")
# 
# #pmse1 <- matrix(0,nrow=nrow(jdf1),ncol=2)
# check1 <- matrix(0,nrow=nrow(jdf1),ncol=length(t_sum))
# check1b <- matrix(0,nrow=nrow(jdf1),ncol=length(t_sum2))
# check1w <- matrix(0,nrow=nrow(jdf1),ncol=length(wsum))
# check1y <- matrix(0,nrow=nrow(jdf1),ncol=length(ysum))
# dist1 <- matrix(0,nrow=nrow(jdf1),ncol=2)
# data1 <- matrix(c(simdata$xee,simdata$xes),ncol=2,byrow=FALSE)
# for(i in 1:nrow(jdf1)){
#   covmat <- matrix(c(jdf1$sigmaxee[i]^2,jdf1$sigmaxee[i]*jdf1$sigmaxes[i]*jdf1$corrx[i],jdf1$sigmaxee[i]*jdf1$sigmaxes[i]*jdf1$corrx[i],jdf1$sigmaxes[i]^2),ncol=2,byrow=TRUE)
#   data <- mvrnorm(length(simdata$xee),c(jdf1$mux1[i],jdf1$mux2[i]),covmat)
#   
#   datawee <- rnorm(length(simdata$xee),data[,1],jdf1$sigmauee[i])
#   datawes <- rnorm(length(simdata$xee),data[,2],jdf1$sigmaues[i])
#   datayee <- rnorm(length(simdata$xee),data[,1],jdf1$sigmaeee[i])
#   datayes <- rnorm(length(simdata$xee),data[,2],jdf1$sigmaees[i])
#   
#   check1[i,1] <- min(data[,1]+data[,2])
#   check1[i,2:7] <- quantile(data[,1]+data[,2],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check1[i,8] <- max(data[,1]+data[,2])
#   check1[i,9] <- min(data[,1])
#   check1[i,10:15] <- quantile(data[,1],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check1[i,16] <- max(data[,1])
#   
#   check1b[i,1] <- sum(data[,1] < 2000)/nrow(data)
#   check1b[i,2] <- sum(data[,2] < 0)/nrow(data)
#   check1b[i,3:4] <- range(2*data[,1] + data[,2])
#   check1b[i,5] <- quantile(data[,1],probs=c(0.75)) - quantile(data[,1],probs=c(0.25))
#   check1b[i,6] <- quantile(data[,1]+data[,2],probs=c(0.75)) - quantile(data[,1]+data[,2],probs=c(0.25))
#   check1b[i,7:8] <- range(data[,2])
#   
#   check1w[i,1] <- min(datawee)
#   check1w[i,2:7] <- quantile(datawee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check1w[i,8] <- max(datawee)
#   check1w[i,9] <- min(datawes)
#   check1w[i,10:15] <- quantile(datawes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check1w[i,16] <- max(datawes)
#   
#   check1y[i,1] <- min(datayee)
#   check1y[i,2:7] <- quantile(datayee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check1y[i,8] <- max(datayee)
#   check1y[i,9] <- min(datayes)
#   check1y[i,10:15] <- quantile(datayes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check1y[i,16] <- max(datayes)
#   
#   dist1[i,1] <- sum(sqrt((data[,1]-simdata$xee)^2+(data[,2]-simdata$xes)^2))
#   dist1[i,2] <- sum(sqrt((data[,1]-data1[,1])^2+(data[,2]-data1[,2])^2))
#   data1 <- data
#   
# }
# 
# pmse1 <- c(mean((simdata$xee - simdata$yee[,1])^2),mean((simdata$xes - simdata$yes[,1])^2))
# m1sum <- colMeans(check1)
# ppred16(check1,t_sum,xnames)
# ppred8(check1b,t_sum2)
# ppred16(check1w,wsum,wnames)
# ppred16(check1y,ysum,wnames)
# 
# check1 <- data.frame(check1)
# qplot(data=check1,x=X2,X3,geom="point") + geom_point(aes(x=t_sum[2],y=t_sum[3]),colour="red",size=5)
# #1,3;1,10;2,10;2,3;
# 
# mdist1 <- melt(data.frame(dist1))
# ggplot(data=mdist1) + geom_histogram(aes(x=value,fill=variable),alpha=0.5,position="")


#--------------------------------------------------
#---------------Model 2----------------------------#
#same as before but allows for linear bias in y
# 
# modelj2 <- "
# model{
# for(i in 1:n){
# for(j in 1:2){
# yee[i,j] ~ dnorm(be0 + be1*x[i,1],taueee)
# yes[i,j] ~ dnorm(bs0 + bs1*x[i,2],tauees)
# wee[i,j] ~ dnorm(x[i,1], tauuee)
# wes[i,j] ~ dnorm(x[i,2], tauues)
# }
# x[i,1:2] ~ dmnorm(mux[],taux[,])
# }
# 
# be0 ~ dnorm(0,1/100000) #dnorm(100,1/10000)
# be1 ~ dnorm(1,1/10000) #dnorm(1,1/10000)
# bs0 ~ dnorm(0,1/100000) #dnorm(0,1/10000)
# bs1 ~ dnorm(1,1/10000) #dnorm(1,1/10000)
# mux[1] ~ dnorm(2200,1/300000)
# mux[2] ~ dnorm(0,1/100000)
# taux[1:2,1:2] ~ dwish(I,d) #d should be dimension + 1
# taueee ~ dgamma(0.1,0.1)#<- (1/sigmaeee)^2 #~ dt(0,1,1)T(0,) half cauchy
# tauees ~ dgamma(0.1,0.1)#<- (1/sigmaees)^2
# tauuee ~ dgamma(0.1,0.1)#<- (1/sigmauee)^2
# tauues ~ dgamma(0.1,0.1)#<- (1/sigmaues)^2
# sigmaeee <- 1/sqrt(taueee)#~ dt(0,1,1)T(0,)
# sigmaees <- 1/sqrt(tauees)#~ dt(0,1,1)T(0,)
# sigmauee <- 1/sqrt(tauuee)#~ dt(0,1,1)T(0,)
# sigmaues <- 1/sqrt(tauues)#~ dt(0,1,1)T(0,)
# sigmax <- inverse(taux[,])
# sigmaxee <- sqrt(sigmax[1,1])
# sigmaxes <- sqrt(sigmax[2,2])
# corrx <- sigmax[1,2]/(sigmaxee*sigmaxes)
# 
# }
# "
# 
# jm2 <- jags.model(textConnection(modelj2),data=simdata,n.adapt=2000,n.chains=3)
# js2 <- coda.samples(jm2,c("be0","be1","bs0","bs1", "mux","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","x"),n.iter=20000)
# #dic2 <- dic.samples(jm2,n.iter=2000)
# #gelman.diag(js2[,c("be0","be1","bs0","bs1","mux[1]","mux[2]", "sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx")],multivariate=FALSE)
# #summary(js2[,c("be0","be1","bs0","bs1", "mux[1]","mux[2]","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx")])
# 
# xeenames <- paste0("x[",1:n,",1]")
# xesnames <- paste0("x[",1:n,",2]")
# 
# 
# m2list <- list(be0=unlist(js2[,"be0"]),be1=unlist(js2[,"be1"]),bs0=unlist(js2[,"bs0"]),
#                bs1=unlist(js2[,"bs1"]),muee=unlist(js2[,"mux[1]"]),mues=unlist(js2[,"mux[2]"]),
#                sigmaeee=unlist(js2[,"sigmaeee"]),sigmaees=unlist(js2[,"sigmaees"]),sigmavee=unlist(js2[,"sigmauee"]),
#                sigmaves=unlist(js2[,"sigmaues"]),sigmaxee=unlist(js2[,"sigmaxee"]),sigmaxes=unlist(js2[,"sigmaxes"]),
#                corrx=unlist(js2[,"corrx"]),latentxee=as.matrix(js2[,xeenames]),
#                latentxes=as.matrix(js2[,xesnames]))
# 
# ppan2 <- pp_jags_lm2(m2list,t_sum,wsum,ysum,simdata$xee,simdata$xes,yeeb,yesb)
# 
#ppred16(ppan2$checkx,t_sum,xnames)
#ppred16(ppan2$checkw,wsum,wnames)
#ppred16(ppan2$checky,ysum,wnames)


# jdf2 <- as.data.frame(as.matrix(js2))
# pmse2 <- matrix(0,nrow=nrow(jdf2),ncol=2)
# names(jdf2)[6:7] <- c("mux1","mux2")
# check2 <- matrix(0,nrow=nrow(jdf2),ncol=length(t_sum))
# check2b <- matrix(0,nrow=nrow(jdf2),ncol=length(t_sum2))
# check2w <- matrix(0,nrow=nrow(jdf2),ncol=length(wsum))
# check2y <- matrix(0,nrow=nrow(jdf2),ncol=length(ysum))
# for(i in 1:nrow(jdf2)){
#   covmat <- matrix(c(jdf2$sigmaxee[i]^2,jdf2$sigmaxee[i]*jdf2$sigmaxes[i]*jdf2$corrx[i],jdf2$sigmaxee[i]*jdf2$sigmaxes[i]*jdf2$corrx[i],jdf2$sigmaxes[i]^2),ncol=2,byrow=TRUE)
#   data <- mvrnorm(length(simdata$xee),c(jdf2$mux1[i],jdf2$mux2[i]),covmat)
#   check2[i,1:7] <- quantile(data[,1]+data[,2],probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95))
#   
#   datawee <- rnorm(length(simdata$xee),data[,1],jdf2$sigmauee[i])
#   datawes <- rnorm(length(simdata$xee),data[,2],jdf2$sigmaues[i])
#   datayee <- rnorm(length(simdata$xee),jdf2$be0[i] + jdf2$be1[i]*data[,1],jdf2$sigmaeee[i])
#   datayes <- rnorm(length(simdata$xee),jdf2$bs0[i] + jdf2$bs1[i]*data[,2],jdf2$sigmaees[i])
#   
#   check2[i,1] <- min(data[,1]+data[,2])
#   check2[i,2:7] <- quantile(data[,1]+data[,2],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check2[i,8] <- max(data[,1]+data[,2])
#   check2[i,9] <- min(data[,1])
#   check2[i,10:15] <- quantile(data[,1],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check2[i,16] <- max(data[,1])
#   
#   check2b[i,1] <- sum(data[,1] < 2000)/nrow(data)
#   check2b[i,2] <- sum(data[,2] < 0)/nrow(data)
#   check2b[i,3:4] <- range(2*data[,1] + data[,2])
#   check2b[i,5] <- quantile(data[,1],probs=c(0.75)) - quantile(data[,1],probs=c(0.25))
#   check2b[i,6] <- quantile(data[,1]+data[,2],probs=c(0.75)) - quantile(data[,1]+data[,2],probs=c(0.25))
#   check2b[i,7:8] <- range(data[,2])
#   
#   check2w[i,1] <- min(datawee)
#   check2w[i,2:7] <- quantile(datawee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check2w[i,8] <- max(datawee)
#   check2w[i,9] <- min(datawes)
#   check2w[i,10:15] <- quantile(datawes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check2w[i,16] <- max(datawes)
#   
#   check2y[i,1] <- min(datayee)
#   check2y[i,2:7] <- quantile(datayee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check2y[i,8] <- max(datayee)
#   check2y[i,9] <- min(datayes)
#   check2y[i,10:15] <- quantile(datayes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check2y[i,16] <- max(datayes)
#   
#   pmse2[i,1] <- mean((jdf2$be0[i] + jdf2$be1[i]*simdata$xee - simdata$yee[,1])^2)
#   pmse2[i,2] <- mean((jdf2$bs0[i] + jdf2$bs1[i]*simdata$xes - simdata$yes[,1])^2)     
# }
# 
# m2sum <- colMeans(check2)
# ppred16(check2,t_sum,xnames)
# ppred8(check2b,t_sum2)
# ppred16(check2w,wsum,wnames)
# ppred16(check2y,ysum,wnames)
# 
# check2 <- data.frame(check2)
# qplot(data=check2,x=X1,X9,geom="point") + geom_point(aes(x=t_sum[1],y=t_sum[9]),colour="red",size=5)

#--------------------------------------------------
#---------------Model 3----------------------------#


modelj3 <- "
model{
for(i in 1:n){
for(j in 1:2){
yee[i,j] ~ dnorm(be0 + be1*x[i,1] + geg*zg[i] + geb*zb[i] + gea*za[i],taueee)
yes[i,j] ~ dnorm(bs0 + bs1*x[i,2] + gig*zg[i] + gib*zb[i] + gia*za[i],tauees)
wee[i,j] ~ dnorm(x[i,1], tauuee)
wes[i,j] ~ dnorm(x[i,2], tauues)
}
x[i,1:2] ~ dmnorm(mux[],taux[,])
}

be0 ~ dnorm(0,1/100000) #dnorm(100,1/10000)
be1 ~ dnorm(1,1/10000) #dnorm(1,1/10000)
bs0 ~ dnorm(0,1/100000) #dnorm(0,1/10000)
bs1 ~ dnorm(1,1/10000) #dnorm(1,1/10000)
mux[1] ~ dnorm(2200,1/300000)
mux[2] ~ dnorm(0,1/100000)
taux[1:2,1:2] ~ dwish(I,d) #d should be dimension + 1
taueee ~ dgamma(0.1,0.1)#<- (1/sigmaeee)^2 #~ dt(0,1,1)T(0,) half cauchy
tauees ~ dgamma(0.1,0.1)#<- (1/sigmaees)^2
tauuee ~ dgamma(0.1,0.1)#<- (1/sigmauee)^2
tauues ~ dgamma(0.1,0.1)#<- (1/sigmaues)^2
sigmaeee <- 1/sqrt(taueee)#~ dt(0,1,1)T(0,)
sigmaees <- 1/sqrt(tauees)#~ dt(0,1,1)T(0,)
sigmauee <- 1/sqrt(tauuee)#~ dt(0,1,1)T(0,)
sigmaues <- 1/sqrt(tauues)#~ dt(0,1,1)T(0,)
sigmax <- inverse(taux[,])
sigmaxee <- sqrt(sigmax[1,1])
sigmaxes <- sqrt(sigmax[2,2])
corrx <- sigmax[1,2]/(sigmaxee*sigmaxes)

geg ~ dnorm(0,1/100000)
geb ~ dnorm(0,1/100000)
gea ~ dnorm(0,1/100000)
gig ~ dnorm(0, 1/100000)
gib ~ dnorm(0,1/100000)
gia ~ dnorm(0,1/100000)
}
"

jm3 <- jags.model(textConnection(modelj3),data=simdata,n.adapt=2000,n.chains=3)
js3 <- coda.samples(jm3,c("be0","be1","bs0","bs1", "mux","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","x","geg","geb","gea","gig","gib","gia"),n.iter=20000)
#dic3 <- dic.samples(jm3,n.iter=2000)
#gelman.diag(js3[,c("be0","be1","bs0","bs1","mux[1]","mux[2]", "sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","geg","geb","gea","gig","gib","gia")],multivariate=FALSE)
#summary(js3[,c("be0","be1","bs0","bs1", "mux[1]","mux[2]","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","geg","geb","gea","gig","gib","gia")])

 xeenames <- paste0("x[",1:n,",1]")
 xesnames <- paste0("x[",1:n,",2]")

m3list <- list(be0=unlist(js3[,"be0"]),be1=unlist(js3[,"be1"]),bs0=unlist(js3[,"bs0"]),
               bs1=unlist(js3[,"bs1"]),muee=unlist(js3[,"mux[1]"]),mues=unlist(js3[,"mux[2]"]),
               sigmaeee=unlist(js3[,"sigmaeee"]),sigmaees=unlist(js3[,"sigmaees"]),sigmavee=unlist(js3[,"sigmauee"]),
               sigmaves=unlist(js3[,"sigmaues"]),sigmaxee=unlist(js3[,"sigmaxee"]),sigmaxes=unlist(js3[,"sigmaxes"]),
               corrx=unlist(js3[,"corrx"]),latentxee=as.matrix(js3[,xeenames]),
               latentxes=as.matrix(js3[,xesnames]),geg=unlist(js3[,"geg"]),
               geb=unlist(js3[,"geb"]),gea=unlist(js3[,"gea"]),gig=unlist(js3[,"gig"]),
               gib=unlist(js3[,"gib"]),gia=unlist(js3[,"gia"]))

ppan3 <- pp_jags_lm(m3list,t_sum,wsum,ysum,simdata$xee,simdata$xes,yeeb,yesb,zg,zb,za)

save.image("//my.files.iastate.edu/Users/dcries/Desktop/jags_mcmc1.RData")

names <- c("$\\beta_{0,ee}$","$\\beta_{1,ee}$","$\\beta_{0,es}$","$\\beta_{1,es}$",
           "$\\mu_x^{EE}$","$\\mu_x^{\\Delta ES}$","$\\sigma_{yee}$","$\\sigma_{yes}$",
           "$\\sigma_{wee}$","$\\sigma_{wes}$","$\\sigma_{xee}$","$\\sigma_{xes}$",
           "$\\rho$","$\\gamma_{1,ee}$","$\\gamma_{2,ee}$","$\\gamma_{3,ee}$",
           "$\\gamma_{1,es}$","$\\gamma_{2,es}$","$\\gamma_{3,es}$")
sumtable <- data.frame(summary(js3[,c("be0","be1","bs0","bs1", "mux[1]","mux[2]","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","geg","geb","gea","gig","gib","gia")])$quantile)
rownames(sumtable) <- names
print(xtable(sumtable), sanitize.text.function=function(x){x})
#ppred16(ppan3$checkx,t_sum,xnames)
#ppred16(ppan3$checkw,wsum,wnames)
#ppred16(ppan3$checky,ysum,wnames)

#------------------------------------------------------------#
#---------------------Model 3 (from sim4.R)------------------#
#instead of single bivariate normal for X, use mixture of 2
# 
# modelj3 <- "
# model{
# #n is number of obs
# #c is number of clusters/mixtures to fit
# #r is number of replicates in observed data
# #2 is for cov(X) is  two dimensional for xee, and xes
# 
# for(i in 1:n){
# for(j in 1:2){
# yee[i,j] ~ dnorm(be0 + be1*x[i,1],1/(sigmaeee^2))
# yes[i,j] ~ dnorm(bs0 + bs1*x[i,2],1/(sigmaees^2))
# wee[i,j] ~ dnorm(x[i,1], 1/(sigmauee^2))
# wes[i,j] ~ dnorm(x[i,2], 1/(sigmaues^2))
# }
# x[i,1:2] ~ dmnorm(mux[zeta[i],1:2],taux[,,zeta[i]])
# zeta[i] ~ dcat(pi[])
# }
# 
# for(k in 1:2){
# mux[k,1:2] ~ dmnorm(mu0[k,1:2],mucov[,,k])
# taux[1:2,1:2,k] ~ dwish(I,d) #d should be dimension + 1
# sigmax[1:2,1:2,k] <- inverse(taux[1:2,1:2,k])
# }
# 
# #  for(i in 1:2){
# #    mux[i,1:2] <- sort(mux0[i,1:2])
# #  }
# 
# sigmaxee <- sqrt(sigmax[1,1,])
# sigmaxes <- sqrt(sigmax[2,2,])
# corrx <- sigmax[1,2,]/(sigmaxee*sigmaxes)
# 
# be0 ~ dnorm(0,1/10000)
# be1 ~ dnorm(1,1/100)
# bs0 ~ dnorm(0,1/10000)
# bs1 ~ dnorm(1,1/100)
# sigmaeee ~ dt(0,1,1)T(0,)
# sigmaees ~ dt(0,1,1)T(0,)
# sigmauee ~ dt(0,1,1)T(0,)
# sigmaues ~ dt(0,1,1)T(0,)
# 
# pi ~ ddirch(a[])
# 
# }
# "
# 
# jm3 <- jags.model(textConnection(modelj3),data=simdata,n.adapt=2000,n.chains=3)
# js3 <- coda.samples(jm3,c("pi", "be0","be1","bs0","bs1", "mux","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","x","zeta"),n.iter=2000)
# dic3 <- dic.samples(jm3,n.iter=2000)
# gelman.diag(js3[,c("be0","be1","bs0","bs1", "sigmaeee","sigmaees","sigmauee","sigmaues","pi[1]","pi[2]")],multivariate=FALSE)
# summary(js3[,c("be0","be1","bs0","bs1","sigmaeee","sigmaees","sigmauee","sigmaues")])
# 
# jdf3 <- as.data.frame(as.matrix(js3))
# pei3 <- jdf3[,21:320] + jdf3[,321:620]
# pmse3 <- sqrt(mean((simdata$xei - colMeans(pei3))^2)) #88.57643
# 
# percentile3 <- as.data.frame(t(apply(pei3,1,quantile,probs=c(0.1,0.25,0.5,0.75,0.9))))
# 
# pmse3 <- matrix(0,nrow=nrow(jdf3),ncol=2)
# names(jdf3)[c(5:10,17:20)] <- c("corrx1","corrx2","mux11","mux21","mux12","mux22","sigmaxee1","sigmaxee2","sigmaxes1","sigmaxes2")
# check3 <- matrix(0,nrow=nrow(jdf3),ncol=length(t_sum))
# data <- matrix(0,nrow=length(simdata$xee),ncol=2)
# for(i in 1:nrow(jdf3)){
#   zeta <- as.numeric(jdf3[i,621:920])
#   
#   covmat1 <- matrix(c(jdf3$sigmaxee1[i]^2,jdf3$sigmaxee1[i]*jdf3$sigmaxes1[i]*jdf3$corrx1[i],jdf3$sigmaxee1[i]*jdf3$sigmaxes1[i]*jdf3$corrx1[i],jdf3$sigmaxes1[i]^2),ncol=2,byrow=TRUE)
#   data1 <- mvrnorm(length(simdata$xee),c(jdf3$mux11[i],jdf3$mux12[i]),covmat1)
#   covmat2 <- matrix(c(jdf3$sigmaxee2[i]^2,jdf3$sigmaxee2[i]*jdf3$sigmaxes2[i]*jdf3$corrx2[i],jdf3$sigmaxee2[i]*jdf3$sigmaxes2[i]*jdf3$corrx2[i],jdf3$sigmaxes2[i]^2),ncol=2,byrow=TRUE)
#   data2 <- mvrnorm(length(simdata$xee),c(jdf3$mux21[i],jdf3$mux22[i]),covmat2)
#   
#   data[which(zeta==1),1:2] <- data1[which(zeta==1),1:2]
#   data[which(zeta==2),1:2] <- data2[which(zeta==2),1:2]
#   
#   check3[i,1:7] <- quantile(data[,1]+data[,2],probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95))
#   check3[i,8:9] <- range(data[,1])
#   check3[i,10] <- sum(data[,1] < 2000)/nrow(data)
#   check3[i,11] <- sum(data[,2] < 0)/nrow(data)
#   check3[i,12:13] <- range(2*data[,1] + data[,2])
#   pmse3[i,1] <- mean((jdf3$be0[i] + jdf3$be1[i]*simdata$xee - simdata$yee[,1])^2)
#   pmse3[i,2] <- mean((jdf3$bs0[i] + jdf3$bs1[i]*simdata$xes - simdata$yes[,1])^2)     
#   
# }
# 
# m3sum <- colMeans(check3)
# ppred(check3,t_sum)
# check3 <- data.frame(check3)
# qplot(data=check3,x=X2,X4,geom="point") + geom_point(aes(x=t_sum[2],y=t_sum[4]),colour="red",size=5)


#--------------------------------------------------------------#
#------------------------Model 4-------------------------------#
#DP
# 
# 
# 
# modelj4 <- "
# model{
# #n is number of obs
# #H is number of clusters/mixtures to fit for DP
# #r is number of replicates in observed data
# #2 is for cov(X) is  two dimensional for xee, and xes
# 
# for(i in 1:n){
# for(j in 1:2){
# yee[i,j] ~ dnorm(be0 + be1*x[i,1],taueee)
# yes[i,j] ~ dnorm(bs0 + bs1*x[i,2],tauees)
# wee[i,j] ~ dnorm(x[i,1],tauuee)
# wes[i,j] ~ dnorm(x[i,2], tauues)
# }
# x[i,1:2] ~ dmnorm(mux[zeta[i],1:2],taux[,,zeta[i]])
# zeta[i] ~ dcat(pi[])
# }
# 
# for(k in 1:H){
# mux[k,1:2] ~ dmnorm(mu0dp,mucovdp)
# taux[1:2,1:2,k] ~ dwish(I,d) #d should be dimension + 1
# sigmax[1:2,1:2,k] <- inverse(taux[1:2,1:2,k])
# }
# 
# for(h in 1:(H-1)){
# V[h] ~ dbeta(1,alpha)
# }
# 
# V[H] <- 1
# pi[1] <- V[1]
# 
# for(h in 2:H){
# pi[h] <- V[h] * (1-V[h-1]) * pi[h-1]/V[h-1]
# }
# 
# #  for(i in 1:c){
# #    mux[i,1:2] <- sort(mux0[i,1:2])
# #  }
# 
# 
# #leave because I will only ever have ee and es
# sigmaxee <- sqrt(sigmax[1,1,])
# sigmaxes <- sqrt(sigmax[2,2,])
# corrx <- sigmax[1,2,]/(sigmaxee*sigmaxes)
# 
# be0 ~ dnorm(0,1/100000)
# be1 ~ dnorm(1,1/10000)
# bs0 ~ dnorm(0,1/100000)
# bs1 ~ dnorm(1,1/10000)
# taueee ~ dgamma(0.1,0.1)#<- (1/sigmaeee)^2 #~ dt(0,1,1)T(0,) half cauchy
# tauees ~ dgamma(0.1,0.1)#<- (1/sigmaees)^2
# tauuee ~ dgamma(0.1,0.1)#<- (1/sigmauee)^2
# tauues ~ dgamma(0.1,0.1)#<- (1/sigmaues)^2
# sigmaeee <- sqrt(taueee)#~ dt(0,1,1)T(0,)
# sigmaees <- sqrt(tauees)#~ dt(0,1,1)T(0,)
# sigmauee <- sqrt(tauuee)#~ dt(0,1,1)T(0,)
# sigmaues <- sqrt(tauues)#~ dt(0,1,1)T(0,)
# 
# #sigmaeee ~ dt(0,1,1)T(0,)
# #sigmaees ~ dt(0,1,1)T(0,)
# #sigmauee ~ dt(0,1,1)T(0,)
# #sigmaues ~ dt(0,1,1)T(0,)
# 
# }
# "
# 
# 
# jm4 <- jags.model(textConnection(modelj4),data=simdata,n.adapt=7000,n.chains=3)
# js4 <- coda.samples(jm4,c("pi", "be0","be1","bs0","bs1", "mux","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","x","zeta"),n.iter=10000)
# dic4 <- dic.samples(jm4,n.iter=2000)
# gelman.diag(js4[,c( "be0","be1","bs0","bs1", "sigmaeee","sigmaees","sigmauee","sigmaues")],multivariate=FALSE)
# summary(js4[,c("be0","be1","bs0","bs1", "sigmaeee","sigmaees","sigmauee","sigmaues")])
# 
# jdf4 <- as.data.frame(as.matrix(js4))
# #order must be be0,be1,bs0,bs1,corrx,mux_1,mux_2,pi,sigmas,sigmaxee,sigmaxes
# pmse4 <- matrix(0,nrow=nrow(jdf4),ncol=2)
# data <- matrix(0,nrow=length(simdata$xee),ncol=2)
# check4 <- matrix(0,nrow=nrow(jdf4),ncol=length(t_sum))
# check4b <- matrix(0,nrow=nrow(jdf4),ncol=length(t_sum2))
# check4w <- matrix(0,nrow=nrow(jdf4),ncol=length(wsum))
# check4y <- matrix(0,nrow=nrow(jdf4),ncol=length(ysum))
# dist4 <- matrix(0,nrow=nrow(jdf4),ncol=2)
# data4 <- matrix(c(simdata$xee,simdata$xes),ncol=2,byrow=FALSE)
# for(i in 1:nrow(jdf4)){
#   zeta <- as.numeric(jdf4[i,729:1028])
#   for(j in 1:20){
#     covmat <- matrix(c(jdf4[i,j+88]^2,jdf4[i,j+88]*jdf4[i,j+108]*jdf4[i,j+4],jdf4[i,j+88]*jdf4[i,j+108]*jdf4[i,j+4],jdf4[i,j+108]^2),ncol=2,byrow=TRUE)
#     datat <- mvrnorm(length(simdata$xee),c(jdf4[i,j+24],jdf4[i,j+44]),covmat)
#     data[which(zeta==j),1:2] <- datat[which(zeta==j),1:2]
#   }
#   
#   datawee <- rnorm(length(simdata$xee),data[,1],jdf4$sigmauee[i])
#   datawes <- rnorm(length(simdata$xee),data[,2],jdf4$sigmaues[i])
#   datayee <- rnorm(length(simdata$xee),jdf4$be0[i] + jdf4$be1[i]*data[,1],jdf4$sigmaeee[i])
#   datayes <- rnorm(length(simdata$xee),jdf4$bs0[i] + jdf4$bs1[i]*data[,2],jdf4$sigmaees[i])
#   
#   check4[i,1] <- min(data[,1]+data[,2])
#   check4[i,2:7] <- quantile(data[,1]+data[,2],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check4[i,8] <- max(data[,1]+data[,2])
#   check4[i,9] <- min(data[,1])
#   check4[i,10:15] <- quantile(data[,1],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check4[i,16] <- max(data[,1])
#   
#   check4b[i,1] <- sum(data[,1] < 2000)/nrow(data)
#   check4b[i,2] <- sum(data[,2] < 0)/nrow(data)
#   check4b[i,3:4] <- range(2*data[,1] + data[,2])
#   check4b[i,5] <- quantile(data[,1],probs=c(0.75)) - quantile(data[,1],probs=c(0.25))
#   check4b[i,6] <- quantile(data[,1]+data[,2],probs=c(0.75)) - quantile(data[,1]+data[,2],probs=c(0.25))
#   check4b[i,7:8] <- range(data[,2])
#   
#   check4w[i,1] <- min(datawee)
#   check4w[i,2:7] <- quantile(datawee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check4w[i,8] <- max(datawee)
#   check4w[i,9] <- min(datawes)
#   check4w[i,10:15] <- quantile(datawes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check4w[i,16] <- max(datawes)
#   
#   check4y[i,1] <- min(datayee)
#   check4y[i,2:7] <- quantile(datayee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check4y[i,8] <- max(datayee)
#   check4y[i,9] <- min(datayes)
#   check4y[i,10:15] <- quantile(datayes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check4y[i,16] <- max(datayes)
#   
#   
#   pmse4[i,1] <- mean((jdf4$be0[i] + jdf4$be1[i]*simdata$xee - simdata$yee[,1])^2)
#   pmse4[i,2] <- mean((jdf4$bs0[i] + jdf4$bs1[i]*simdata$xes - simdata$yes[,1])^2)     
#   
#   dist4[i,1] <- sum(sqrt((data[,1]-simdata$xee)^2+(data[,2]-simdata$xes)^2))
#   dist4[i,2] <- sum(sqrt((data[,1]-data4[,1])^2+(data[,2]-data4[,2])^2))
#   data4 <- data
#   
#   
# }
# ppred16(check4,t_sum,xnames)
# ppred8(check4b,t_sum2)
# ppred16(check4w,wsum,wnames)
# ppred16(check4y,ysum,wnames)
# 
# m4sum <- colMeans(check4)
# check4 <- data.frame(check4)
# qplot(data=check4,x=X2,X4,geom="point") + geom_point(aes(x=t_sum[2],y=t_sum[4]),colour="red",size=5)
# 
# mdist4 <- melt(data.frame(dist4))
# ggplot(data=mdist4) + geom_histogram(aes(x=value,fill=variable),alpha=0.5,position="")
# 
# 
# loc_check <- seq(1300,3300,length=10)
# density_conv(jdf4[,c(25:44,65:84,89:108)],loc_check[1])
# 
# loc_check <- seq(-1000,1000,length=10)
# density_conv(jdf4[,c(45:64,65:84,109:128)],loc_check[1])
# 
# table(apply(jdf4[,729:ncol(jdf4)],1,function(x){return(length(unique(x)))}))
# 

#-------------------------------------------------#
#---------------------Model 5---------------------#
# 
# modelj5 <- "
# model{
# #n is number of obs
# #H is number of clusters/mixtures to fit for DP
# #r is number of replicates in observed data
# #2 is for cov(X) is  two dimensional for xee, and xes
# 
# for(i in 1:n){
# for(j in 1:2){
# yee[i,j] ~ dnorm(be0 + be1*x[i,1] + geg*zg[i] + geb*zb[i] + gea*za[i],taueee)
# yes[i,j] ~ dnorm(bs0 + bs1*x[i,2] + gig*zg[i] + gib*zb[i] + gia*za[i],tauees)
# wee[i,j] ~ dnorm(x[i,1], tauuee)
# wes[i,j] ~ dnorm(x[i,2], tauues)
# }
# x[i,1:2] ~ dmnorm(mux[zeta[i],1:2],taux[,,zeta[i]])
# zeta[i] ~ dcat(pi[])
# }
# 
# for(k in 1:H){
# mux[k,1:2] ~ dmnorm(mu0dp,mucovdp)
# taux[1:2,1:2,k] ~ dwish(I,d) #d should be dimension + 1
# sigmax[1:2,1:2,k] <- inverse(taux[1:2,1:2,k])
# }
# 
# for(h in 1:(H-1)){
# V[h] ~ dbeta(1,alpha)
# }
# 
# 
# V[H] <- 1
# pi[1] <- V[1]
# 
# for(h in 2:H){
# pi[h] <- V[h] * (1-V[h-1]) * pi[h-1]/V[h-1]
# }
# 
# #  for(i in 1:c){
# #    mux[i,1:2] <- sort(mux0[i,1:2])
# #  }
# 
# 
# #leave because I will only ever have ee and es
# sigmaxee <- sqrt(sigmax[1,1,])
# sigmaxes <- sqrt(sigmax[2,2,])
# corrx <- sigmax[1,2,]/(sigmaxee*sigmaxes)
# 
# be0 ~ dnorm(0,1/100000)
# be1 ~ dnorm(1,1/10000)
# bs0 ~ dnorm(0,1/100000)
# bs1 ~ dnorm(1,1/10000)
# #sigmaeee ~ dt(0,1,1)T(0,)
# #sigmaees ~ dt(0,1,1)T(0,)
# #sigmauee ~ dt(0,1,1)T(0,)
# #sigmaues ~ dt(0,1,1)T(0,)
# taueee ~ dgamma(0.1,0.1)#<- (1/sigmaeee)^2 #~ dt(0,1,1)T(0,) half cauchy
# tauees ~ dgamma(0.1,0.1)#<- (1/sigmaees)^2
# tauuee ~ dgamma(0.1,0.1)#<- (1/sigmauee)^2
# tauues ~ dgamma(0.1,0.1)#<- (1/sigmaues)^2
# sigmaeee <- 1/sqrt(taueee)#~ dt(0,1,1)T(0,)
# sigmaees <- 1/sqrt(tauees)#~ dt(0,1,1)T(0,)
# sigmauee <- 1/sqrt(tauuee)#~ dt(0,1,1)T(0,)
# sigmaues <- 1/sqrt(tauues)#~ dt(0,1,1)T(0,)
# 
# 
# geg ~ dnorm(0,1/100000)
# geb ~ dnorm(0,1/100000)
# gea ~ dnorm(0,1/100000)
# gig ~ dnorm(0, 1/100000)
# gib ~ dnorm(0,1/100000)
# gia ~ dnorm(0,1/100000)
# 
# }
# "
# 
# jm5 <- jags.model(textConnection(modelj5),data=simdata,n.adapt=2000,n.chains=3)
# js5 <- coda.samples(jm5,c( "be0","be1","bs0","bs1", "mux","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx",  "x","geg","geb","gea","gig","gib","gia","zeta","pi"),n.iter=10000)
# #dic5 <- dic.samples(jm5,n.iter=2000)
# #gelman.diag(js5[,c("be0","be1","bs0","bs1","sigmaeee","sigmaees","sigmauee","sigmaues","geg","geb","gea","gig","gib","gia")],multivariate=FALSE)
# #summary(js5[,c("be0","be1","bs0","bs1", "sigmaeee","sigmaees","sigmauee","sigmaues","geg","geb","gea","gig","gib","gia")])
# 
# 
# mueenames <- paste0("mux[",1:simdata$H,",1]")
# muesnames <- paste0("mux[",1:simdata$H,",2]")
# sxeenames <- paste0("sigmaxee[",1:simdata$H,"]")
# sxesnames <- paste0("sigmaxes[",1:simdata$H,"]")
# cornames <- paste0("corrx[",1:simdata$H,"]")
# zetanames <- paste0("zeta[",1:n,"]")
# pinames <- paste0("pi[",1:simdata$H,"]")
# 
# 
# pisum = as.matrix(js5[,pinames])
# 
# m5list <- list(be0=unlist(js5[,"be0"]),be1=unlist(js5[,"be1"]),bs0=unlist(js5[,"bs0"]),
#                bs1=unlist(js5[,"bs1"]),muee=as.matrix(js5[,mueenames]),mues=as.matrix(js5[,muesnames]),
#                sigmaeee=unlist(js5[,"sigmaeee"]),sigmaees=unlist(js5[,"sigmaees"]),sigmavee=unlist(js5[,"sigmauee"]),
#                sigmaves=unlist(js5[,"sigmaues"]),sigmaxee=as.matrix(js5[,sxeenames]),sigmaxes=as.matrix(js5[,sxesnames]),
#                corrx=as.matrix(js5[,cornames]),latentxee=as.matrix(js5[,xeenames]),
#                latentxes=as.matrix(js5[,xesnames]),geg=unlist(js5[,"geg"]),
#                geb=unlist(js5[,"geb"]),gea=unlist(js5[,"gea"]),gig=unlist(js5[,"gig"]),
#                gib=unlist(js5[,"gib"]),gia=unlist(js5[,"gia"]),zeta=as.matrix(js5[,zetanames]))
# 
# ppan5 <- pp_jags_dp(m5list,t_sum,wsum,ysum,simdata$xee,simdata$xes,yeeb,yesb,zg,zb,za)
# 
# ppred16(ppan5$checkx,t_sum,xnames)
# ppred16(ppan5$checkw,wsum,wnames)
# ppred16(ppan5$checky,ysum,wnames)

save.image("//my.files.iastate.edu/Users/dcries/Desktop/jags_mcmc3.RData")


# jdf5 <- as.data.frame(as.matrix(js5))
# data <- matrix(0,nrow=length(simdata$xee),ncol=2)
# pmse5 <- matrix(0,nrow=nrow(jdf5),ncol=2)
# check5 <- matrix(0,nrow=nrow(jdf5),ncol=length(t_sum))
# check5b <- matrix(0,nrow=nrow(jdf5),ncol=length(t_sum2))
# check5w <- matrix(0,nrow=nrow(jdf5),ncol=length(wsum))
# check5y <- matrix(0,nrow=nrow(jdf5),ncol=length(ysum))
# for(i in 1:nrow(jdf5)){
#   zeta <- as.numeric(jdf5[i,735:1034])
#   for(j in 1:20){
#     covmat <- matrix(c(jdf5[i,j+94]^2,jdf5[i,j+94]*jdf5[i,j+114]*jdf5[i,j+4],jdf5[i,j+94]*jdf5[i,j+114]*jdf5[i,j+4],jdf5[i,j+114]^2),ncol=2,byrow=TRUE)
#     datat <- mvrnorm(length(simdata$xee),c(jdf5[i,j+30],jdf5[i,j+50]),covmat)
#     data[which(zeta==j),1:2] <- datat[which(zeta==j),1:2]
#   }
#   
#   datawee <- rnorm(length(simdata$xee),data[,1],jdf5$sigmauee[i])
#   datawes <- rnorm(length(simdata$xee),data[,2],jdf5$sigmaues[i])
#   datayee <- rnorm(length(simdata$xee),jdf5$be0[i] + jdf5$be1[i]*data[,1] + zg*jdf5$geg[i] + zb*jdf5$geb[i] + za*jdf5$gea[i],jdf5$sigmaeee[i])
#   datayes <- rnorm(length(simdata$xee),jdf5$bs0[i] + jdf5$bs1[i]*data[,2] + zg*jdf5$gig[i] + zb*jdf5$gib[i] + za*jdf5$gia[i],jdf5$sigmaees[i])
#   
#   
#   check5[i,1] <- min(data[,1]+data[,2])
#   check5[i,2:7] <- quantile(data[,1]+data[,2],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check5[i,8] <- max(data[,1]+data[,2])
#   check5[i,9] <- min(data[,1])
#   check5[i,10:15] <- quantile(data[,1],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check5[i,16] <- max(data[,1])
#   
#   check5b[i,1] <- sum(data[,1] < 2000)/nrow(data)
#   check5b[i,2] <- sum(data[,2] < 0)/nrow(data)
#   check5b[i,3:4] <- range(2*data[,1] + data[,2])
#   check5b[i,5] <- quantile(data[,1],probs=c(0.75)) - quantile(data[,1],probs=c(0.25))
#   check5b[i,6] <- quantile(data[,1]+data[,2],probs=c(0.75)) - quantile(data[,1]+data[,2],probs=c(0.25))
#   check5b[i,7:8] <- range(data[,2])
#   
#   check5w[i,1] <- min(datawee)
#   check5w[i,2:7] <- quantile(datawee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check5w[i,8] <- max(datawee)
#   check5w[i,9] <- min(datawes)
#   check5w[i,10:15] <- quantile(datawes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check5w[i,16] <- max(datawes)
#   
#   check5y[i,1] <- min(datayee)
#   check5y[i,2:7] <- quantile(datayee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check5y[i,8] <- max(datayee)
#   check5y[i,9] <- min(datayes)
#   check5y[i,10:15] <- quantile(datayes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check5y[i,16] <- max(datayes)
#   
#   pmse5[i,1] <- mean((jdf5$be0[i] + jdf5$be1[i]*simdata$xee + zg*jdf5$geg[i] + zb*jdf5$geb[i] + za*jdf5$gea[i] - simdata$yee[,1])^2)
#   pmse5[i,2] <- mean((jdf5$bs0[i] + jdf5$bs1[i]*simdata$xes + zg*jdf5$gig[i] + zb*jdf5$gib[i] + za*jdf5$gia[i] - simdata$yes[,1])^2)     
#   
# }
# ppred16(check5,t_sum,xnames)
# ppred8(check5b,t_sum2)
# ppred16(check5w,wsum,wnames)
# ppred16(check5y,ysum,wnames)
# 
# check5 <- data.frame(check5)
# qplot(data=check5,x=X2,X4,geom="point") + geom_point(aes(x=t_sum[2],y=t_sum[4]),colour="red",size=5)
# 
# 
# 
# loc_check <- seq(1300,3300,length=10)
# density_conv(jdf5[,c(31:50,71:90,95:114)],loc_check[1])
# 
# loc_check <- seq(-1000,1000,length=10)
# density_conv(jdf5[,c(51:70,71:90,115:134)],loc_check[1])
# 
# table(apply(jdf5[,735:ncol(jdf5)],1,function(x){return(length(unique(x)))}))

#---------------------------------------------------------------------------#
#--------------------------Model 6------------------------------------------#
#DP, B-spline + Z covariates for Y
# 
# modelj6 <- "
# model{
# #n is number of obs
# #H is number of clusters/mixtures to fit for DP
# #r is number of replicates in observed data
# #2 is for cov(X) is  two dimensional for xee, and xes
# 
# for(i in 1:n){
# meanee[i] <- inprod(B[i,],beta[]) 
# meanes[i] <- inprod(B2[i,],beta2[]) 
# for(j in 1:2){
# yee[i,j] ~ dnorm(meanee[i] + geg*zg[i] + geb*zb[i] + gea*za[i],1/(sigmaeee^2))
# yes[i,j] ~ dnorm(meanes[i]+ gig*zg[i] + gib*zb[i] + gia*za[i],1/(sigmaees^2))
# wee[i,j] ~ dnorm(x[i,1], 1/(sigmauee^2))
# wes[i,j] ~ dnorm(x[i,2], 1/(sigmaues^2))
# }
# x[i,1:2] ~ dmnorm(mux[zeta[i],1:2],taux[,,zeta[i]])
# zeta[i] ~ dcat(pi[])
# 
# rescxm[i,1] <- (x[i,1]-knots[1])/(knots[2]-knots[1])
# rescxm[i,2] <- (x[i,1]-knots[2])/(knots[3]-knots[2])
# rescxm[i,K] <- (x[i,1]-knots[K-1])/(knots[K]-knots[K-1])
# B[i,1] <- ifelse((x[i,1]>=knots[1]) && (x[i,1]<=knots[2]),(1/2)*(1-rescxm[i,1])^2,0) 
# check[i] <- ifelse((x[i,1]>=knots[2])&&(x[i,1]<=knots[3]),1,0)
# val1[i] <- ifelse((check[i]==0) && (B[i,1]!=0), -(rescxm[i,1]^2)+rescxm[i,1]+1/2, 0)
# B[i,2] <- ifelse(check[i]==1,(1/2)*(1-rescxm[i,2])^2,val1[i])
# for(j in 3:(K-1)){
# rescxm[i,j] <- (x[i,1]-knots[j])/(knots[j+1]-knots[j])
# check1m[i,j] <- ifelse((x[i,1]>=knots[j])&&(x[i,1]<=knots[j+1]),1,0)
# check2m[i,j] <- ifelse((x[i,1]>=knots[j-1])&&(x[i,1]<=knots[j]),1,0)
# check3m[i,j] <- ifelse((x[i,1]>=knots[j-2])&&(x[i,1]<=knots[j-1]),1,0)
# val1m[i,j] <- ifelse(check3m[i,j]==1, (rescxm[i,j-2]^2)/2,0)
# val2m[i,j] <- ifelse(check2m[i,j]==1, -(rescxm[i,j-1]^2)+rescxm[i,j-1]+1/2,val1m[i,j])
# B[i,j] <- ifelse(check1m[i,j]==1, (1/2)*(1-rescxm[i,j])^2,val2m[i,j])
# }
# check2[i] <- ifelse((x[i,1]>=knots[K-2])&&(x[i,1]<=knots[K-1]),1,0)
# check3[i] <- ifelse((x[i,1]>=knots[K-1])&&(x[i,1]<=knots[K]),1,0)
# check4[i] <- ifelse(check2[i]==1 || check3[i]==1,1,0)
# val2[i] <- ifelse(check2[i]==1, (rescxm[i,K-2]^2)/2,-(rescxm[i,K-1]^2)+rescxm[i,K-1]+1/2)
# B[i,K] <- ifelse(check4[i]==0,0,val2[i])
# B[i,K+1] <- ifelse((x[i,1]>=knots[K-1])&&(x[i,1]<=knots[K]),(rescxm[i,K-1]^2)/2,0)
# 
# 
# rescxm2[i,1] <- (x[i,2]-knots2[1])/(knots2[2]-knots2[1])
# rescxm2[i,2] <- (x[i,2]-knots2[2])/(knots2[3]-knots2[2])
# rescxm2[i,K] <- (x[i,2]-knots2[K-1])/(knots2[K]-knots2[K-1])
# B2[i,1] <- ifelse((x[i,2]>=knots2[1]) && (x[i,2]<=knots2[2]),(1/2)*(1-rescxm2[i,1])^2,0) 
# checka[i] <- ifelse((x[i,2]>=knots2[2])&&(x[i,2]<=knots2[3]),1,0)
# val12[i] <- ifelse((checka[i]==0) && (B2[i,1]!=0), -(rescxm2[i,1]^2)+rescxm2[i,1]+1/2, 0)
# B2[i,2] <- ifelse(checka[i]==1,(1/2)*(1-rescxm2[i,2])^2,val12[i])
# for(j in 3:(K-1)){
# rescxm2[i,j] <- (x[i,2]-knots2[j])/(knots2[j+1]-knots2[j])
# check1m2[i,j] <- ifelse((x[i,2]>=knots2[j])&&(x[i,2]<=knots2[j+1]),1,0)
# check2m2[i,j] <- ifelse((x[i,2]>=knots2[j-1])&&(x[i,2]<=knots2[j]),1,0)
# check3m2[i,j] <- ifelse((x[i,2]>=knots2[j-2])&&(x[i,2]<=knots2[j-1]),1,0)
# val1m2[i,j] <- ifelse(check3m2[i,j]==1, (rescxm2[i,j-2]^2)/2,0)
# val2m2[i,j] <- ifelse(check2m2[i,j]==1, -(rescxm2[i,j-1]^2)+rescxm2[i,j-1]+1/2,val1m2[i,j])
# B2[i,j] <- ifelse(check1m2[i,j]==1, (1/2)*(1-rescxm2[i,j])^2,val2m2[i,j])
# }
# check22[i] <- ifelse((x[i,2]>=knots2[K-2])&&(x[i,2]<=knots2[K-1]),1,0)
# check32[i] <- ifelse((x[i,2]>=knots2[K-1])&&(x[i,2]<=knots2[K]),1,0)
# check42[i] <- ifelse(check22[i]==1 || check32[i]==1,1,0)
# val22[i] <- ifelse(check22[i]==1, (rescxm2[i,K-2]^2)/2,-(rescxm2[i,K-1]^2)+rescxm2[i,K-1]+1/2)
# B2[i,K] <- ifelse(check42[i]==0,0,val22[i])
# B2[i,K+1] <- ifelse((x[i,2]>=knots2[K-1])&&(x[i,2]<=knots2[K]),(rescxm2[i,K-1]^2)/2,0)
# 
# }
# 
# for(i in 1:(K+1)){
# #beta[i] ~ dnorm(0,1/1000)
# beta[i] ~ dunif(-1000000,1000000)
# beta2[i] ~ dunif(-1000000,1000000)
# }
# 
# for(k in 1:H){
# mux[k,1:2] ~ dmnorm(mu0dp,mucovdp)
# taux[1:2,1:2,k] ~ dwish(I,d) #d should be dimension + 1
# sigmax[1:2,1:2,k] <- inverse(taux[1:2,1:2,k])
# }
# 
# for(h in 1:(H-1)){
# V[h] ~ dbeta(1,alpha)
# }
# 
# V[H] <- 1
# pi[1] <- V[1]
# 
# for(h in 2:H){
# pi[h] <- V[h] * (1-V[h-1]) * pi[h-1]/V[h-1]
# }
# 
# 
# #leave because I will only ever have ee and es
# sigmaxee <- sqrt(sigmax[1,1,])
# sigmaxes <- sqrt(sigmax[2,2,])
# corrx <- sigmax[1,2,]/(sigmaxee*sigmaxes)
# 
# sigmaeee ~ dt(0,1,1)T(0,)
# sigmaees ~ dt(0,1,1)T(0,)
# sigmauee ~ dt(0,1,1)T(0,)
# sigmaues ~ dt(0,1,1)T(0,)
# 
# geg ~ dnorm(0,1/100000)
# geb ~ dnorm(0,1/100000)
# gea ~ dnorm(0,1/100000)
# gig ~ dnorm(0,1/100000)
# gib ~ dnorm(0,1/100000)
# gia ~ dnorm(0,1/100000)
# 
# }
# "
# 
# jm6 <- jags.model(textConnection(modelj6),data=simdata,n.adapt=500000,n.chains=5) #18000 is enough for adaptation
# js6 <- coda.samples(jm6,c( "mux","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx",  "x","geg","geb","gea","gig","gib","gia","meanee","meanes","zeta","beta","beta2","pi"),n.iter=5000)
# dic6 <- dic.samples(jm6,n.iter=20000)
# gelman.diag(js6[,c("sigmaeee","sigmaees","sigmauee","sigmaues","geg","geb","gea","gig","gib","gia")],multivariate=FALSE)
# summary(js6[,c( "sigmaeee","sigmaees","sigmauee","sigmaues","geg","geb","gea","gig","gib","gia")])
# 
# jdf6 <- as.data.frame(as.matrix(js6))
# 
# data <- matrix(0,nrow=length(simdata$xee),ncol=2)
# pmse6 <- matrix(0,nrow=nrow(jdf6),ncol=2)
# check6 <- matrix(0,nrow=nrow(jdf6),ncol=length(t_sum))
# check6b <- matrix(0,nrow=nrow(jdf6),ncol=length(t_sum2))
# basis1 <- B.basis(simdata$xee,simdata$knots)
# basis2 <- B.basis(simdata$xes,simdata$knots2)
# for(i in 1:nrow(jdf6)){
#   zeta <- as.numeric(jdf6[i,1347:1646])
#   for(j in 1:20){
#     covmat <- matrix(c(jdf6[i,j+706]^2,jdf6[i,j+706]*jdf6[i,j+16]*jdf6[i,j+726],jdf6[i,j+706]*jdf6[i,j+16]*jdf6[i,j+726],jdf6[i,j+726]^2),ncol=2,byrow=TRUE)
#     datat <- mvrnorm(length(simdata$xee),c(jdf6[i,j+642],jdf6[i,j+662]),covmat)
#     data[which(zeta==j),1:2] <- datat[which(zeta==j),1:2]
#   }
#   
#   check6[i,1] <- min(data[,1]+data[,2])
#   check6[i,2:7] <- quantile(data[,1]+data[,2],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check6[i,8] <- max(data[,1]+data[,2])
#   check6[i,9] <- min(data[,1])
#   check6[i,10:15] <- quantile(data[,1],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check6[i,16] <- max(data[,1])
#   
#   check6b[i,1] <- sum(data[,1] < 2000)/nrow(data)
#   check6b[i,2] <- sum(data[,2] < 0)/nrow(data)
#   check6b[i,3:4] <- range(2*data[,1] + data[,2])
#   check6b[i,5] <- quantile(data[,1],probs=c(0.75)) - quantile(data[,1],probs=c(0.25))
#   check6b[i,6] <- quantile(data[,1]+data[,2],probs=c(0.75)) - quantile(data[,1]+data[,2],probs=c(0.25))
#   check6b[i,7:8] <- range(data[,2])
#   
#   #basis1 <- B.basis(jdf6[i,727:1026],simdata$knots)
#   #basis2 <- B.basis(jdf6[i,1027:1326],simdata$knots2)
#   #pmse6[i,1] <- mean((basis1%*%as.numeric(jdf6[i,1:8]) + zg*jdf6$geg[i] + zb*jdf6$geb[i] + za*jdf6$gea[i] - simdata$yee[,1])^2)
#   #pmse6[i,2] <- mean((basis2%*%as.numeric(jdf6[i,9:16]) + zg*jdf6$gig[i] + zb*jdf6$gib[i] + za*jdf6$gia[i] - simdata$yes[,1])^2)     
#   pmse6[i,1] <- mean((B.basis(as.numeric(jdf6[i,747:1046]),simdata$knots)%*%as.numeric(jdf6[i,1:8]) + zg*jdf6$geg[i] + zb*jdf6$geb[i] + za*jdf6$gea[i] - simdata$yee[,1])^2)
#   pmse6[i,2] <- mean((B.basis(as.numeric(jdf6[i,1047:1346]),simdata$knots2)%*%as.numeric(jdf6[i,9:16]) + zg*jdf6$gig[i] + zb*jdf6$gib[i] + za*jdf6$gia[i] - simdata$yes[,1])^2)     
#   #print(i)
# }
# 
# ppred16(check6,t_sum)
# ppred8(check6b,t_sum2)
# 
# ggplot() + geom_line(aes(x=as.vector(unlist(jdf6[5590,747:1046])),y=as.vector(unlist(jdf6[5590,43:342]))),colour="red") + geom_point(aes(x=simdata$wee[,1],y=simdata$yee[,1])) + geom_line(aes(x=simdata$xee,y=eecurve(simdata$xee)),colour="yellow")#+ geom_line(aes(x=X,y=beta0*sin(beta1*X)),colour="yellow")
# ggplot() + geom_line(aes(x=as.vector(unlist(jdf6[940,1047:1346])),y=as.vector(unlist(jdf6[940,343:642]))),colour="red") + geom_point(aes(x=simdata$wes[,1],y=simdata$yes[,1]))  + geom_line(aes(x=simdata$xes,y=escurve(simdata$xes)),colour="yellow")#+ geom_line(aes(x=X,y=beta0*sin(beta1*X)),colour="yellow")
# 
# 
# loc_check <- seq(1300,3300,length=10)
# density_conv(jdf6[,c(643:662,683:702,707:726)],loc_check[1])
# 
# loc_check <- seq(-1000,1000,length=10)
# density_conv(jdf6[,c(663:682,683:702,727:746)],loc_check[1])
# 
# table(apply(jdf6[,1347:ncol(jdf6)],1,function(x){return(length(unique(x)))}))

#----------------------------------------------------------------#
#-------------------Model 7--------------------------------------#
#b-spline for Y + linear in known Z, X is assumed only single MVN
# 
# modelj7 <- "
# #var B[n,(K+1),niter];
# model
# {
#   
#   for(i in 1:n){
#   mean[i] <- inprod(B[i,],beta[]) #+ mu
#   mean2[i] <- inprod(B2[i,],beta2[]) #+ mu
#   x[i,1:2] ~ dmnorm(mux[],taux[,])
#   
#   for(j in 1:2){
#   yee[i,j] ~ dnorm(mean[i]+geg*zg[i]+gea*za[i]+geb*zb[i], 1/sqrt(sigmaeee))
#   yes[i,j] ~ dnorm(mean2[i]+gig*zg[i]+gia*za[i]+gib*zb[i], 1/sqrt(sigmaees))
#   wee[i,j] ~ dnorm(x[i,1],1/sqrt(sigmauee))
#   wes[i,j] ~ dnorm(x[i,2],1/sqrt(sigmaues))
#   
#   }
#   
#   rescxm[i,1] <- (x[i,1]-knots[1])/(knots[2]-knots[1])
#   rescxm[i,2] <- (x[i,1]-knots[2])/(knots[3]-knots[2])
#   rescxm[i,K] <- (x[i,1]-knots[K-1])/(knots[K]-knots[K-1])
#   B[i,1] <- ifelse((x[i,1]>=knots[1]) && (x[i,1]<=knots[2]),(1/2)*(1-rescxm[i,1])^2,0) 
#   check[i] <- ifelse((x[i,1]>=knots[2])&&(x[i,1]<=knots[3]),1,0)
#   val1[i] <- ifelse((check[i]==0) && (B[i,1]!=0), -(rescxm[i,1]^2)+rescxm[i,1]+1/2, 0)
#   B[i,2] <- ifelse(check[i]==1,(1/2)*(1-rescxm[i,2])^2,val1[i])
#   for(j in 3:(K-1)){
#   rescxm[i,j] <- (x[i,1]-knots[j])/(knots[j+1]-knots[j])
#   check1m[i,j] <- ifelse((x[i,1]>=knots[j])&&(x[i,1]<=knots[j+1]),1,0)
#   check2m[i,j] <- ifelse((x[i,1]>=knots[j-1])&&(x[i,1]<=knots[j]),1,0)
#   check3m[i,j] <- ifelse((x[i,1]>=knots[j-2])&&(x[i,1]<=knots[j-1]),1,0)
#   val1m[i,j] <- ifelse(check3m[i,j]==1, (rescxm[i,j-2]^2)/2,0)
#   val2m[i,j] <- ifelse(check2m[i,j]==1, -(rescxm[i,j-1]^2)+rescxm[i,j-1]+1/2,val1m[i,j])
#   B[i,j] <- ifelse(check1m[i,j]==1, (1/2)*(1-rescxm[i,j])^2,val2m[i,j])
#   }
#   check2[i] <- ifelse((x[i,1]>=knots[K-2])&&(x[i,1]<=knots[K-1]),1,0)
#   check3[i] <- ifelse((x[i,1]>=knots[K-1])&&(x[i,1]<=knots[K]),1,0)
#   check4[i] <- ifelse(check2[i]==1 || check3[i]==1,1,0)
#   val2[i] <- ifelse(check2[i]==1, (rescxm[i,K-2]^2)/2,-(rescxm[i,K-1]^2)+rescxm[i,K-1]+1/2)
#   B[i,K] <- ifelse(check4[i]==0,0,val2[i])
#   B[i,K+1] <- ifelse((x[i,1]>=knots[K-1])&&(x[i,1]<=knots[K]),(rescxm[i,K-1]^2)/2,0)
#   
#   
#   rescxm2[i,1] <- (x[i,2]-knots2[1])/(knots2[2]-knots2[1])
#   rescxm2[i,2] <- (x[i,2]-knots2[2])/(knots2[3]-knots2[2])
#   rescxm2[i,K] <- (x[i,2]-knots2[K-1])/(knots2[K]-knots2[K-1])
#   B2[i,1] <- ifelse((x[i,2]>=knots2[1]) && (x[i,2]<=knots2[2]),(1/2)*(1-rescxm2[i,1])^2,0) 
#   checka[i] <- ifelse((x[i,2]>=knots2[2])&&(x[i,2]<=knots2[3]),1,0)
#   val12[i] <- ifelse((checka[i]==0) && (B2[i,1]!=0), -(rescxm2[i,1]^2)+rescxm2[i,1]+1/2, 0)
#   B2[i,2] <- ifelse(checka[i]==1,(1/2)*(1-rescxm2[i,2])^2,val12[i])
#   for(j in 3:(K-1)){
#   rescxm2[i,j] <- (x[i,2]-knots2[j])/(knots2[j+1]-knots2[j])
#   check1m2[i,j] <- ifelse((x[i,2]>=knots2[j])&&(x[i,2]<=knots2[j+1]),1,0)
#   check2m2[i,j] <- ifelse((x[i,2]>=knots2[j-1])&&(x[i,2]<=knots2[j]),1,0)
#   check3m2[i,j] <- ifelse((x[i,2]>=knots2[j-2])&&(x[i,2]<=knots2[j-1]),1,0)
#   val1m2[i,j] <- ifelse(check3m2[i,j]==1, (rescxm2[i,j-2]^2)/2,0)
#   val2m2[i,j] <- ifelse(check2m2[i,j]==1, -(rescxm2[i,j-1]^2)+rescxm2[i,j-1]+1/2,val1m2[i,j])
#   B2[i,j] <- ifelse(check1m2[i,j]==1, (1/2)*(1-rescxm2[i,j])^2,val2m2[i,j])
#   }
#   check22[i] <- ifelse((x[i,2]>=knots2[K-2])&&(x[i,2]<=knots2[K-1]),1,0)
#   check32[i] <- ifelse((x[i,2]>=knots2[K-1])&&(x[i,2]<=knots2[K]),1,0)
#   check42[i] <- ifelse(check22[i]==1 || check32[i]==1,1,0)
#   val22[i] <- ifelse(check22[i]==1, (rescxm2[i,K-2]^2)/2,-(rescxm2[i,K-1]^2)+rescxm2[i,K-1]+1/2)
#   B2[i,K] <- ifelse(check42[i]==0,0,val22[i])
#   B2[i,K+1] <- ifelse((x[i,2]>=knots2[K-1])&&(x[i,2]<=knots2[K]),(rescxm2[i,K-1]^2)/2,0)
#   
#   }
#   
#   for(i in 1:(K+1)){
#   #beta[i] ~ dnorm(0,1/10000000000)
#   beta[i] ~ dunif(-1000000,1000000)
#   beta2[i] ~ dunif(-1000000,1000000)
#   }
#   
#   sigmaeee ~ dt(0,1,1)T(0,)
#   sigmaees ~ dt(0,1,1)T(0,)
#   sigmauee ~ dt(0,1,1)T(0,)
#   sigmaues ~ dt(0,1,1)T(0,)
#   
#   geg ~ dnorm(0,1/100000)
#   geb ~ dnorm(0,1/100000)
#   gea ~ dnorm(0,1/100000)
#   gig ~ dnorm(0, 1/100000)
#   gib ~ dnorm(0,1/100000)
#   gia ~ dnorm(0,1/100000)
#   
#   mux[1] ~ dnorm(2000,1/10000)
#   mux[2] ~ dnorm(0,1/10000)
#   
#   taux[1:2,1:2] ~ dwish(I,3) #d should be dimension + 1
#   sigmax <- inverse(taux[,])
#   sigmaxee <- sqrt(sigmax[1,1])
#   sigmaxes <- sqrt(sigmax[2,2])
#   corrx <- sigmax[1,2]/(sigmaxee*sigmaxes)
#   
#   
# }
# "
# 
# 
# jm7 <- jags.model(textConnection(modelj7),data=simdata,n.adapt=20000,n.chains=3)
# js7 <- coda.samples(jm7,c( "mux","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx",  "x","mean","mean2","beta","beta2","geg","geb","gea","gig","gib","gia"),n.iter=20000)
# dic7 <- dic.samples(jm7,n.iter=20000)
# gelman.diag(js7[,c("sigmaeee","sigmaees","sigmauee","sigmaues","geg","geb","gea","gig","gib","gia")],multivariate=FALSE)
# summary(js7[,c( "sigmaeee","sigmaees","sigmauee","sigmaues","geg","geb","gea","gig","gib","gia")])
# 
# jdf7 <- as.data.frame(as.matrix(js7))
# names(jdf7)[624:625] <- c("mux1","mux2")
# 
# pmse7 <- matrix(0,nrow=nrow(jdf7),ncol=2)
# check7 <- matrix(0,nrow=nrow(jdf7),ncol=length(t_sum))
# check7b <- matrix(0,nrow=nrow(jdf7),ncol=length(t_sum2))
# basis1 <- B.basis(simdata$xee,simdata$knots)
# basis2 <- B.basis(simdata$xes,simdata$knots2)
# for(i in 1:nrow(jdf7)){
#   covmat <- matrix(c(jdf7$sigmaxee[i]^2,jdf7$sigmaxee[i]*jdf7$sigmaxes[i]*jdf7$corrx[i],jdf7$sigmaxee[i]*jdf7$sigmaxes[i]*jdf7$corrx[i],jdf7$sigmaxes[i]^2),ncol=2,byrow=TRUE)
#   data <- mvrnorm(length(simdata$xee),c(jdf7$mux1[i],jdf7$mux2[i]),covmat)
#   
#   check7[i,1] <- min(data[,1]+data[,2])
#   check7[i,2:7] <- quantile(data[,1]+data[,2],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check7[i,8] <- max(data[,1]+data[,2])
#   check7[i,9] <- min(data[,1])
#   check7[i,10:15] <- quantile(data[,1],probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
#   check7[i,16] <- max(data[,1])
#   
#   check7b[i,1] <- sum(data[,1] < 2000)/nrow(data)
#   check7b[i,2] <- sum(data[,2] < 0)/nrow(data)
#   check7b[i,3:4] <- range(2*data[,1] + data[,2])
#   check7b[i,5] <- quantile(data[,1],probs=c(0.75)) - quantile(data[,1],probs=c(0.25))
#   check7b[i,6] <- quantile(data[,1]+data[,2],probs=c(0.75)) - quantile(data[,1]+data[,2],probs=c(0.25))
#   check7b[i,7:8] <- range(data[,2])
#   
#   #basis1 <- B.basis(jdf7[i,727:1027],simdata$knots)
#   #basis2 <- B.basis(jdf7[i,1027:1327],simdata$knots2)
#   #pmse7[i,1] <- mean((basis1%*%as.numeric(jdf7[i,1:8]) + zg*jdf7$geg[i] + zb*jdf7$geb[i] + za*jdf7$gea[i] - simdata$yee[,1])^2)
#   #pmse7[i,2] <- mean((basis2%*%as.numeric(jdf7[i,9:16]) + zg*jdf7$gig[i] + zb*jdf7$gib[i] + za*jdf7$gia[i] - simdata$yes[,1])^2)     
#   pmse7[i,1] <- mean((B.basis(as.numeric(jdf7[i,632:931]),simdata$knots)%*%as.numeric(jdf7[i,1:8]) + zg*jdf7$geg[i] + zb*jdf7$geb[i] + za*jdf7$gea[i] - simdata$yee[,1])^2)
#   pmse7[i,2] <- mean((B.basis(as.numeric(jdf7[i,932:1231]),simdata$knots2)%*%as.numeric(jdf7[i,9:16]) + zg*jdf7$gig[i] + zb*jdf7$gib[i] + za*jdf7$gia[i] - simdata$yes[,1])^2)     
#   
# }
# 
# ppred16(check7,t_sum)
# ppred8(check7b,t_sum2)
# 
# 
# ggplot() + geom_line(aes(x=as.vector(unlist(jdf7[11400,632:931])),y=as.vector(unlist(jdf7[11400,24:323]))),colour="red") + geom_point(aes(x=simdata$wee[,1],y=simdata$yee[,1])) + geom_line(aes(x=simdata$xee,y=eecurve(simdata$xee)),colour="yellow")#+ geom_line(aes(x=X,y=beta0*sin(beta1*X)),colour="yellow")
# ggplot() + geom_line(aes(x=as.vector(unlist(jdf7[11400,932:1231])),y=as.vector(unlist(jdf7[11400,324:623]))),colour="red") + geom_point(aes(x=simdata$wes[,1],y=simdata$yes[,1]))  + geom_line(aes(x=simdata$xes,y=escurve(simdata$xes)),colour="yellow")#+ geom_line(aes(x=X,y=beta0*sin(beta1*X)),colour="yellow")


#----------------------------------------------------------------#
#-------------------Model 8--------------------------------------#
#b-spline for Y + linear in known Z, X is assumed only single MVN

# modelj8 <- "
# #var B[n,(K+1),niter];
# model
# {
#   
#   for(i in 1:n){
#   mean[i] <- inprod(B[i,],beta[]) #+ mu
#   mean2[i] <- inprod(B2[i,],beta2[]) #+ mu
#   x[i,1:2] ~ dmnorm(mux[,i],taux[,])
#   mux[1,i] ~ dnorm(ae + aeg*zg[i] + aeb*zb[i] + aea*za[i],1/10000)
#   mux[2,i] ~ dnorm(as + asg*zg[i] + asb*zb[i] + asa*za[i],1/10000)
#   
#   for(j in 1:2){
#   yee[i,j] ~ dnorm(mean[i]+geg*zg[i]+gea*za[i]+geb*zb[i], 1/sqrt(sigmaeee))
#   yes[i,j] ~ dnorm(mean2[i]+gig*zg[i]+gia*za[i]+gib*zb[i], 1/sqrt(sigmaees))
#   wee[i,j] ~ dnorm(x[i,1],1/sqrt(sigmauee))
#   wes[i,j] ~ dnorm(x[i,2],1/sqrt(sigmaues))
#   
#   }
#   
#   rescxm[i,1] <- (x[i,1]-knots[1])/(knots[2]-knots[1])
#   rescxm[i,2] <- (x[i,1]-knots[2])/(knots[3]-knots[2])
#   rescxm[i,K] <- (x[i,1]-knots[K-1])/(knots[K]-knots[K-1])
#   B[i,1] <- ifelse((x[i,1]>=knots[1]) && (x[i,1]<=knots[2]),(1/2)*(1-rescxm[i,1])^2,0) 
#   check[i] <- ifelse((x[i,1]>=knots[2])&&(x[i,1]<=knots[3]),1,0)
#   val1[i] <- ifelse((check[i]==0) && (B[i,1]!=0), -(rescxm[i,1]^2)+rescxm[i,1]+1/2, 0)
#   B[i,2] <- ifelse(check[i]==1,(1/2)*(1-rescxm[i,2])^2,val1[i])
#   for(j in 3:(K-1)){
#   rescxm[i,j] <- (x[i,1]-knots[j])/(knots[j+1]-knots[j])
#   check1m[i,j] <- ifelse((x[i,1]>=knots[j])&&(x[i,1]<=knots[j+1]),1,0)
#   check2m[i,j] <- ifelse((x[i,1]>=knots[j-1])&&(x[i,1]<=knots[j]),1,0)
#   check3m[i,j] <- ifelse((x[i,1]>=knots[j-2])&&(x[i,1]<=knots[j-1]),1,0)
#   val1m[i,j] <- ifelse(check3m[i,j]==1, (rescxm[i,j-2]^2)/2,0)
#   val2m[i,j] <- ifelse(check2m[i,j]==1, -(rescxm[i,j-1]^2)+rescxm[i,j-1]+1/2,val1m[i,j])
#   B[i,j] <- ifelse(check1m[i,j]==1, (1/2)*(1-rescxm[i,j])^2,val2m[i,j])
#   }
#   check2[i] <- ifelse((x[i,1]>=knots[K-2])&&(x[i,1]<=knots[K-1]),1,0)
#   check3[i] <- ifelse((x[i,1]>=knots[K-1])&&(x[i,1]<=knots[K]),1,0)
#   check4[i] <- ifelse(check2[i]==1 || check3[i]==1,1,0)
#   val2[i] <- ifelse(check2[i]==1, (rescxm[i,K-2]^2)/2,-(rescxm[i,K-1]^2)+rescxm[i,K-1]+1/2)
#   B[i,K] <- ifelse(check4[i]==0,0,val2[i])
#   B[i,K+1] <- ifelse((x[i,1]>=knots[K-1])&&(x[i,1]<=knots[K]),(rescxm[i,K-1]^2)/2,0)
#   
#   
#   rescxm2[i,1] <- (x[i,2]-knots2[1])/(knots2[2]-knots2[1])
#   rescxm2[i,2] <- (x[i,2]-knots2[2])/(knots2[3]-knots2[2])
#   rescxm2[i,K] <- (x[i,2]-knots2[K-1])/(knots2[K]-knots2[K-1])
#   B2[i,1] <- ifelse((x[i,2]>=knots2[1]) && (x[i,2]<=knots2[2]),(1/2)*(1-rescxm2[i,1])^2,0) 
#   checka[i] <- ifelse((x[i,2]>=knots2[2])&&(x[i,2]<=knots2[3]),1,0)
#   val12[i] <- ifelse((checka[i]==0) && (B2[i,1]!=0), -(rescxm2[i,1]^2)+rescxm2[i,1]+1/2, 0)
#   B2[i,2] <- ifelse(checka[i]==1,(1/2)*(1-rescxm2[i,2])^2,val12[i])
#   for(j in 3:(K-1)){
#   rescxm2[i,j] <- (x[i,2]-knots2[j])/(knots2[j+1]-knots2[j])
#   check1m2[i,j] <- ifelse((x[i,2]>=knots2[j])&&(x[i,2]<=knots2[j+1]),1,0)
#   check2m2[i,j] <- ifelse((x[i,2]>=knots2[j-1])&&(x[i,2]<=knots2[j]),1,0)
#   check3m2[i,j] <- ifelse((x[i,2]>=knots2[j-2])&&(x[i,2]<=knots2[j-1]),1,0)
#   val1m2[i,j] <- ifelse(check3m2[i,j]==1, (rescxm2[i,j-2]^2)/2,0)
#   val2m2[i,j] <- ifelse(check2m2[i,j]==1, -(rescxm2[i,j-1]^2)+rescxm2[i,j-1]+1/2,val1m2[i,j])
#   B2[i,j] <- ifelse(check1m2[i,j]==1, (1/2)*(1-rescxm2[i,j])^2,val2m2[i,j])
#   }
#   check22[i] <- ifelse((x[i,2]>=knots2[K-2])&&(x[i,2]<=knots2[K-1]),1,0)
#   check32[i] <- ifelse((x[i,2]>=knots2[K-1])&&(x[i,2]<=knots2[K]),1,0)
#   check42[i] <- ifelse(check22[i]==1 || check32[i]==1,1,0)
#   val22[i] <- ifelse(check22[i]==1, (rescxm2[i,K-2]^2)/2,-(rescxm2[i,K-1]^2)+rescxm2[i,K-1]+1/2)
#   B2[i,K] <- ifelse(check42[i]==0,0,val22[i])
#   B2[i,K+1] <- ifelse((x[i,2]>=knots2[K-1])&&(x[i,2]<=knots2[K]),(rescxm2[i,K-1]^2)/2,0)
#   
#   }
#   
#   for(i in 1:(K+1)){
#   #beta[i] ~ dnorm(0,1/10000000000)
#   beta[i] ~ dunif(-1000000,1000000)
#   beta2[i] ~ dunif(-1000000,1000000)
#   }
#   
#   sigmaeee ~ dt(0,1,1)T(0,)
#   sigmaees ~ dt(0,1,1)T(0,)
#   sigmauee ~ dt(0,1,1)T(0,)
#   sigmaues ~ dt(0,1,1)T(0,)
#   
#   geg ~ dnorm(0,1/100000)
#   geb ~ dnorm(0,1/100000)
#   gea ~ dnorm(0,1/100000)
#   gig ~ dnorm(0, 1/100000)
#   gib ~ dnorm(0,1/100000)
#   gia ~ dnorm(0,1/100000)
#   
#   ae ~ dnorm(0,1/100000)
#   aeg ~ dnorm(0,1/100000)
#   aeb ~ dnorm(0,1/100000)
#   aea ~ dnorm(0,1/100000)
#   asg ~ dnorm(0, 1/100000)
#   asb ~ dnorm(0,1/100000)
#   asa ~ dnorm(0,1/100000)
#   as ~ dnorm(0,1/100000)
#   
#   taux[1:2,1:2] ~ dwish(I,3) #d should be dimension + 1
#   sigmax <- inverse(taux[,])
#   sigmaxee <- sqrt(sigmax[1,1])
#   sigmaxes <- sqrt(sigmax[2,2])
#   corrx <- sigmax[1,2]/(sigmaxee*sigmaxes)
#   
#   
# }
# "
# 
# 
# jm8 <- jags.model(textConnection(modelj8),data=simdata,n.adapt=2000,n.chains=3)
# js8 <- coda.samples(jm8,c( "mux","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx",  "x","mean","mean2","beta","beta2","geg","geb","gea","gig","gib","gia","ae","aeg","aeb","aea","as","asg","asb","asa"),n.iter=2000)
# dic8 <- dic.samples(jm8,n.iter=2000)
# gelman.diag(js8[,c("sigmaeee","sigmaees","sigmauee","sigmaues","geg","geb","gea","gig","gib","gia","ae","aeg","aeb","aea","as","asg","asb","asa")],multivariate=FALSE)
# summary(js8[,c( "sigmaeee","sigmaees","sigmauee","sigmaues","geg","geb","gea","gig","gib","gia","ae","aeg","aeb","aea","as","asg","asb","asa")])
# 
# jdf8 <- as.data.frame(as.matrix(js8))
# names(jdf8)[624:625] <- c("mux1","mux2")
# 
# pmse8 <- matrix(0,nrow=nrow(jdf8),ncol=2)
# check8 <- matrix(0,nrow=nrow(jdf8),ncol=length(t_sum))
# basis1 <- B.basis(simdata$xee,simdata$knots)
# basis2 <- B.basis(simdata$xes,simdata$knots2)
# for(i in 1:nrow(jdf8)){
#   covmat <- matrix(c(jdf8$sigmaxee[i]^2,jdf8$sigmaxee[i]*jdf8$sigmaxes[i]*jdf8$corrx[i],jdf8$sigmaxee[i]*jdf8$sigmaxes[i]*jdf8$corrx[i],jdf8$sigmaxes[i]^2),ncol=2,byrow=TRUE)
#   data <- mvrnorm(length(simdata$xee),c(jdf8$mux1[i],jdf8$mux2[i]),covmat)
#   check8[i,1:7] <- quantile(data[,1]+data[,2],probs=c(0.05,0.1,0.25,0.5,0.85,0.9,0.95))
#   check8[i,8:9] <- range(data[,1])
#   check8[i,10] <- sum(data[,1] < 2000)/nrow(data)
#   check8[i,11] <- sum(data[,2] < 0)/nrow(data)
#   check8[i,12:13] <- range(2*data[,1] + data[,2])
#   #basis1 <- B.basis(jdf8[i,828:1028],simdata$knots)
#   #basis2 <- B.basis(jdf8[i,1028:1328],simdata$knots2)
#   #pmse8[i,1] <- mean((basis1%*%as.numeric(jdf8[i,1:8]) + zg*jdf8$geg[i] + zb*jdf8$geb[i] + za*jdf8$gea[i] - simdata$yee[,1])^2)
#   #pmse8[i,2] <- mean((basis2%*%as.numeric(jdf8[i,9:16]) + zg*jdf8$gig[i] + zb*jdf8$gib[i] + za*jdf8$gia[i] - simdata$yes[,1])^2)     
#   pmse8[i,1] <- mean((B.basis(as.numeric(jdf8[i,632:931]),simdata$knots)%*%as.numeric(jdf8[i,1:8]) + zg*jdf8$geg[i] + zb*jdf8$geb[i] + za*jdf8$gea[i] - simdata$yee[,1])^2)
#   pmse8[i,2] <- mean((B.basis(as.numeric(jdf8[i,932:1231]),simdata$knots2)%*%as.numeric(jdf8[i,9:16]) + zg*jdf8$gig[i] + zb*jdf8$gib[i] + za*jdf8$gia[i] - simdata$yes[,1])^2)     
#   
# }
# 
# ppred(check8,t_sum)
# 

#ggplot() + geom_line(aes(x=as.vector(unlist(jdf7[400,609:908])),y=as.vector(unlist(jdf7[400,1:300]))),colour="red") + geom_point(aes(x=simdata$wee[,1],y=simdata$yee[,1])) + geom_line(aes(x=simdata$xee,y=eecurve(simdata$xee)),colour="yellow")#+ geom_line(aes(x=X,y=beta0*sin(beta1*X)),colour="yellow")
#ggplot() + geom_line(aes(x=as.vector(unlist(jdf7[400,909:1208])),y=as.vector(unlist(jdf7[400,301:600]))),colour="red") + geom_point(aes(x=simdata$wes[,1],y=simdata$yes[,1]))  + geom_line(aes(x=simdata$xes,y=escurve(simdata$xes)),colour="yellow")#+ geom_line(aes(x=X,y=beta0*sin(beta1*X)),colour="yellow")
