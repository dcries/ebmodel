library(rstan)
library(rjags)
library(MASS)
library(ggplot2)
library(gridExtra)
library(mgcv) 
library(reshape)
library(mvtnorm)
library(fGarch)

Rcpp::sourceCpp('C:/Users/dcries/github/ebmodel/ppred_analysis.cpp')
source('C:\\Users\\dcries\\github\\ebmodel\\base_fcn.R')
setwd('/home/dcries')
Rcpp::sourceCpp('ebmodel/ppred_analysis.cpp')
source('ebmodel/base_fcn.R')

set.seed(936)
#number of simulated data sets
nsim <- 50
#number of posterior predictive iterations
m <- 1000
#burn in
nadapt <- 2000
#iterations
niter <- 10000
#n chains
nchain <- 1
#pick which individuals interested in
ind <- c(72,188,259)
nind <- length(ind)

pmse0 <- matrix(0,nrow=2,ncol=nsim)
pmse1 <- matrix(0,nrow=2,ncol=nsim)

model0 <- matrix(0,nrow=12,ncol=nsim)
model1 <- matrix(0,nrow=19,ncol=nsim)

indcheck <- matrix(0,nrow=nsim,ncol=4*length(ind))
indnames <- c(paste0("x[",ind[1],",1]"),paste0("x[",ind[1],",2]"),
              paste0("x[",ind[2],",1]"),paste0("x[",ind[2],",2]"),
              paste0("x[",ind[3],",1]"),paste0("x[",ind[3],",2]"))

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


params <- c(100,50,300,14,-7,-200,8,-5)
#names for X diagnostics
wnames <- c("EE Min","EE 5%tile","EE 10%tile","EE 25%tile","EE 75%tile","EE 90%tile","EE 95%tile","EE Max",
            "ES Min","ES 5%tile","ES 10%tile","ES 25%tile","ES 75%tile","ES 90%tile","ES 95%tile","ES Max")


for(i in 1:nsim){
  simdata <- generate_data5(params,dist=1)
  simdata$d <- 3 # for inv-wish prior
  simdata$I <- diag(2) #for inv-wish prior
  simdata$n <- length(simdata$xee)
  
  mu0 <- matrix(c(2800,0,2400,0,2000,0),ncol=2,byrow=TRUE)
  mucov <- array(0,dim=c(2,2,3))
  diag(mucov[,,1]) <- c(1/(200^2),1/(200^2))
  diag(mucov[,,2]) <- c(1/(200^2),1/(200^2))
  diag(mucov[,,3]) <- c(1/(200^2),1/(200^2))
  
  simdata$mu0 <- mu0
  simdata$mucov <- mucov

  
  #quantiles of xei, range of xee, P(xee < 2000), P(xes < 0), ie burning more than gaining,
  #t_sum <- c(quantile(simdata$xei,probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)),range(simdata$xee),sum(simdata$xee < 2000)/length(simdata$xee),sum(simdata$xes < 0)/length(simdata$xes),range(simdata$xee+simdata$xei))
  t_sum <- c(min(simdata$xee),quantile(simdata$xee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$xee),min(simdata$xes),quantile(simdata$xes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$xes))
  #discrepancy measures for W measurements
  wsum <- c(min(simdata$wee),quantile(simdata$wee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$wee),min(simdata$wes),quantile(simdata$wes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$wes))
  #discrepancy meaures for Y measurements
  ysum <- c(min(simdata$yee),quantile(simdata$yee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$yee),min(simdata$yes),quantile(simdata$yes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$yes))
  
  
  
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
  
  
  
  naivedata <- list(yee=simdata$yee,yes=simdata$yes,x=matrix(c(weeb,wesb),ncol=2,byrow=FALSE),n=length(simdata$xee),I=diag(2),d=3,zg=zg,za=za,zb=zb)
  
  jm0 <- jags.model(textConnection(modelj0),data=naivedata,n.adapt=nadapt,n.chains=nchain,quiet=TRUE)
  js0 <- coda.samples(jm0,c("mux","sigmaeee","sigmaees","be0","be1","bs0","bs1","geg","geb","gea","gig","gib","gia"),n.iter=niter,quiet=TRUE)
  jdf0 <- as.data.frame(as.matrix(js0))
  
  pp_indicies <- sample(1:(nrow(jdf0)/3),m)
  #jdf0 <- jdf0[pp_indicies,]
  
  check0y <- matrix(0,nrow=m,ncol=length(wsum))
  dist0 <- matrix(0,nrow=m,ncol=2)
  data0 <- matrix(c(simdata$xee,simdata$xes),ncol=2,byrow=FALSE)
  pmse0a <- matrix(0,ncol=2,nrow=m)
  for(j in 1:m){
    datayee <- rnorm(length(simdata$xee),jdf0$be0[j]+jdf0$be1[j]*weeb,jdf0$sigmaeee[j])
    datayes <- rnorm(length(simdata$xee),jdf0$bs0[j]+jdf0$bs1[j]*wesb,jdf0$sigmaees[j])
    
    check0y[j,1] <- min(datayee)
    check0y[j,2:7] <- quantile(datayee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
    check0y[j,8] <- max(datayee)
    check0y[j,9] <- min(datayes)
    check0y[j,10:15] <- quantile(datayes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95))
    check0y[j,16] <- max(datayes)
    
    pmse0a[j,1] <- mean((yeeb - jdf0$be0[j]-jdf0$be1[j]*weeb-jdf0$geg[j]*zg-jdf0$geb[j]*zb-jdf0$gea[j]*za)^2)
    pmse0a[j,2] <- mean((yesb - jdf0$bs0[j]-jdf0$bs1[j]*wesb-jdf0$gig[j]*zg-jdf0$gib[j]*zb-jdf0$gia[j]*za)^2)
  }
  
  #sumtable0 <- data.frame(summary(js0[,c("be0","be1","bs0","bs1", "sigmaeee","sigmaees","geg","geb","gea","gig","gib","gia")])$quantile)
  
  
  jm3 <- jags.model(textConnection(modelj3),data=simdata,n.adapt=nadapt,n.chains=nchain,quiet=TRUE)
  js3 <- coda.samples(jm3,c("be0","be1","bs0","bs1", "mux","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","x","geg","geb","gea","gig","gib","gia"),n.iter=niter,quiet=TRUE)
  
  sumtable <- data.frame(summary(js3[,c("be0","be1","bs0","bs1", "mux[1]","mux[2]","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","geg","geb","gea","gig","gib","gia")])$quantile)
  
  xeenames <- paste0("x[",1:n,",1]")
  xesnames <- paste0("x[",1:n,",2]")
  
  #js3 <- js3[pp_indicies,]
  m3list <- list(be0=unlist(js3[,"be0"])[pp_indicies],be1=unlist(js3[,"be1"])[pp_indicies],bs0=unlist(js3[,"bs0"])[pp_indicies],
                 bs1=unlist(js3[,"bs1"])[pp_indicies],muee=unlist(js3[,"mux[1]"])[pp_indicies],mues=unlist(js3[,"mux[2]"])[pp_indicies],
                 sigmaeee=unlist(js3[,"sigmaeee"])[pp_indicies],sigmaees=unlist(js3[,"sigmaees"])[pp_indicies],sigmavee=unlist(js3[,"sigmauee"])[pp_indicies],
                 sigmaves=unlist(js3[,"sigmaues"])[pp_indicies],sigmaxee=unlist(js3[,"sigmaxee"])[pp_indicies],sigmaxes=unlist(js3[,"sigmaxes"])[pp_indicies],
                 corrx=unlist(js3[,"corrx"])[pp_indicies],latentxee=as.matrix(js3[,xeenames])[pp_indicies,],
                 latentxes=as.matrix(js3[,xesnames])[pp_indicies,],geg=unlist(js3[,"geg"])[pp_indicies],
                 geb=unlist(js3[,"geb"])[pp_indicies],gea=unlist(js3[,"gea"])[pp_indicies],gig=unlist(js3[,"gig"])[pp_indicies],
                 gib=unlist(js3[,"gib"])[pp_indicies],gia=unlist(js3[,"gia"])[pp_indicies])
  
  ppan3 <- pp_jags_lm(m3list,t_sum,wsum,ysum,simdata$xee,simdata$xes,yeeb,yesb,zg,zb,za)
  
  
  pmse0[,i] <- colMeans(pmse0a)
  pmse1[,i] <- colMeans(ppan3$pmse)
  model0[,i] <- summary(js0)$statistics[,1]
  model1[,i] <- summary(js3[,c("be0","be1","bs0","bs1", "mux[1]","mux[2]","sigmaeee","sigmaees","sigmauee","sigmaues","sigmaxee","sigmaxes","corrx","geg","geb","gea","gig","gib","gia")])$statistics[,1]

  indcheck[i,1] <- simdata$xee[ind[1]]
  indcheck[i,2] <- mean(unlist(js3[,indnames[1]]))
  indcheck[i,3] <- simdata$xes[ind[1]]
  indcheck[i,4] <- mean(unlist(js3[,indnames[2]]))
  
  indcheck[i,5] <- simdata$xee[ind[2]]
  indcheck[i,6] <- mean(unlist(js3[,indnames[3]]))
  indcheck[i,7] <- simdata$xes[ind[2]]
  indcheck[i,8] <- mean(unlist(js3[,indnames[4]]))
  
  indcheck[i,9] <- simdata$xee[ind[3]]
  indcheck[i,10] <- mean(unlist(js3[,indnames[5]]))
  indcheck[i,11] <- simdata$xes[ind[3]]
  indcheck[i,12] <- mean(unlist(js3[,indnames[6]]))

}

pmse0 <- t(pmse0)
pmse1 <- t(pmse1)
model0 <- t(model0)
model1 <- t(model1)

out <- list(pmse0=pmse0,pmse1=pmse1,model0=model0,model1=model1)
save(out,file="\\\\my.files.iastate.edu\\Users\\dcries\\Desktop\\simple_2_1.RData")

