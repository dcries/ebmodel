library(ggplot2)
library(coda)
library(reshape)
library(splines)
library(gridExtra)
library(mvtnorm)
library(MASS)
library(fGarch)


source('//my.files.iastate.edu/Users/dcries/Desktop/research/rprograms/base_fcn.R', echo=FALSE)#Rcpp::sourceCpp('//my.files.iastate.edu/Users/dcries/Desktop/research/rprograms/trial_and_simulated/bivar_fullmcmc2.cpp')
load("//my.files.iastate.edu/Users/dcries/Desktop/bvn_mcmc2x.RData")
params <- c(100,50,300,14,-7,-200,8,-5)
simdata <- generate_data5(params,dist=2)

yee <- simdata$yee
yes <- simdata$yes
yeeb <- rowMeans(yee)
yesb <- rowMeans(yes)
xee <- simdata$xee
xes <- simdata$xes
wee <- simdata$wee
wes <- simdata$wes
n <- length(yeeb)
# 
zg <- simdata$zg #<- rbinom(n,1,0.5) #gender indicator
zb <- simdata$zb#<- rnorm(n,27,5) #bmi
za <- simdata$za#<- runif(n,20,40) #age
Z= cbind(zg,zb,za)

chain1 <- out2x$chain1
chain2 <- out2x$chain2
chain3 <- out2x$chain3


latentxee <- mcmc.list(mcmc(chain1$latentxee),mcmc(chain2$latentxee),mcmc(chain3$latentxee))
latentxes <- mcmc.list(mcmc(chain1$latentxes),mcmc(chain2$latentxes),mcmc(chain3$latentxes))
muee <- mcmc.list(mcmc(chain1$muee),mcmc(chain2$muee),mcmc(chain3$muee))
mues <- mcmc.list(mcmc(chain1$mues),mcmc(chain2$mues),mcmc(chain3$mues))
sigma2xee <- mcmc.list(mcmc(chain1$sigma2xee),mcmc(chain2$sigma2xee),mcmc(chain3$sigma2xee))
sigma2xes <- mcmc.list(mcmc(chain1$sigma2xes),mcmc(chain2$sigma2xes),mcmc(chain3$sigma2xes))
corrx <- mcmc.list(mcmc(chain1$corrx),mcmc(chain2$corrx),mcmc(chain3$corrx))
meanfcnee <- mcmc.list(mcmc(chain1$meanfcnee),mcmc(chain2$meanfcnee),mcmc(chain3$meanfcnee))
meanfcnes <- mcmc.list(mcmc(chain1$meanfcnes),mcmc(chain2$meanfcnes),mcmc(chain3$meanfcnes))
sigma2eee <- mcmc.list(mcmc(chain1$sigma2eee),mcmc(chain2$sigma2eee),mcmc(chain3$sigma2eee))
sigma2ees <- mcmc.list(mcmc(chain1$sigma2ees),mcmc(chain2$sigma2ees),mcmc(chain3$sigma2ees))
sigma2vee <- mcmc.list(mcmc(chain1$sigma2vee),mcmc(chain2$sigma2vee),mcmc(chain3$sigma2vee))
sigma2ves <- mcmc.list(mcmc(chain1$sigma2ves),mcmc(chain2$sigma2ves),mcmc(chain3$sigma2ves))
ree <- mcmc.list(mcmc(chain1$ree),mcmc(chain2$ree),mcmc(chain3$ree))
res <- mcmc.list(mcmc(chain1$res),mcmc(chain2$res),mcmc(chain3$res))
betaee <- mcmc.list(mcmc(chain1$betaee),mcmc(chain2$betaee),mcmc(chain3$betaee))
betaes <- mcmc.list(mcmc(chain1$betaes),mcmc(chain2$betaes),mcmc(chain3$betaes))
gammaee <- mcmc.list(mcmc(chain1$gammaee),mcmc(chain2$gammaee),mcmc(chain3$gammaee))
gammaes <- mcmc.list(mcmc(chain1$gammaes),mcmc(chain2$gammaes),mcmc(chain3$gammaes))
kee <- mcmc.list(mcmc(apply(chain1$ree,1,function(x) length(unique(x)))-1),mcmc(apply(chain2$ree,1,function(x) length(unique(x)))-1),mcmc(apply(chain3$ree,1,function(x) length(unique(x)))-1))
kes <- mcmc.list(mcmc(apply(chain1$res,1,function(x) length(unique(x)))-1),mcmc(apply(chain2$res,1,function(x) length(unique(x)))-1),mcmc(apply(chain3$res,1,function(x) length(unique(x)))-1))



plot(xee,apply(chain1$latentxes,2,function(x){length(unique(x))/length(x)}))
plot(xee,apply(chain2$latentxes,2,function(x){length(unique(x))/length(x)}))
plot(xee,apply(chain3$latentxes,2,function(x){length(unique(x))/length(x)}))
chain1$acceptance_rate

gelman.diag(muee)
gelman.diag(mues)
gelman.diag(sigma2xee)
gelman.diag(sigma2xes)
gelman.diag(corrx)
gelman.diag(sigma2eee)
gelman.diag(sigma2ees)
gelman.diag(sigma2vee)
gelman.diag(sigma2ves)
gelman.diag(kee)
gelman.diag(kes)
gelman.diag(betaee)
gelman.diag(betaes)

plot(betaee)
pairs(chain1$betaee)
pairs(data.frame(chain1$sigma2eee,chain1$sigma2ees,chain1$sigma2vee,chain1$sigma2ves))

truth <- c(min(simdata$xee),quantile(simdata$xee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$xee),min(simdata$xes),quantile(simdata$xes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$xes))
truth2 <- c(sum(simdata$xee < 2000)/length(simdata$xee),sum(simdata$xes < 0)/length(simdata$xes),range(simdata$xee+simdata$xei),quantile(simdata$xee,probs=c(0.75))-quantile(simdata$xee,probs=c(0.25)),quantile(simdata$xei,probs=c(0.75))-quantile(simdata$xei,probs=c(0.25)),range(simdata$xes))
wsum <- c(min(simdata$wee),quantile(simdata$wee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$wee),min(simdata$wes),quantile(simdata$wes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$wes))
ysum <- c(min(simdata$yee),quantile(simdata$yee,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$yee),min(simdata$yes),quantile(simdata$yes,probs=c(0.05,0.1,0.25,0.75,0.9,0.95)),max(simdata$yes))


ppan=pp_bvn(chain1,truth,wsum,ysum,xee,xes,yeeb,yesb,Z)

out2x <- list(chain1=chain1,chain2=chain2,chain3=chain3,ppan=ppan)
save(out2x,file="//my.files.iastate.edu/Users/dcries/Desktop/bvn_mcmc2x.RData")
#save.image("//my.files.iastate.edu/Users/dcries/Desktop/full_mcmc1x.RData")


# xnames <- c("EI Min","EI 5%tile","EI 10%tile","EI 25%tile","EI 75%tile","EI 90%tile","EI 95%tile","EI Max",
#             "EE Min","EE 5%tile","EE 10%tile","EE 25%tile","EE 75%tile","EE 90%tile","EE 95%tile","EE Max")

wnames <- c("EE Min","EE 5%tile","EE 10%tile","EE 25%tile","EE 75%tile","EE 90%tile","EE 95%tile","EE Max",
            "ES Min","ES 5%tile","ES 10%tile","ES 25%tile","ES 75%tile","ES 90%tile","ES 95%tile","ES Max")

ppred16(ppan$checkx,truth,wnames)
ppred16(ppan$checkw,wsum,wnames)
ppred16(ppan$checky,ysum,wnames)


summat <- matrix(0,nrow=15,ncol=5)
summat[1,] <- summary(muee)$quantiles
summat[2,] <- summary(mues)$quantiles
summat[3,] <- quantile(sqrt(unlist(sigma2eee)),probs=c(0.025,0.25,0.5,0.75,0.975))
summat[4,] <- quantile(sqrt(unlist(sigma2ees)),probs=c(0.025,0.25,0.5,0.75,0.975))
summat[5,] <- quantile(sqrt(unlist(sigma2vee)),probs=c(0.025,0.25,0.5,0.75,0.975))
summat[6,] <- quantile(sqrt(unlist(sigma2ves)),probs=c(0.025,0.25,0.5,0.75,0.975))
summat[7,] <- quantile(sqrt(unlist(sigma2xee)),probs=c(0.025,0.25,0.5,0.75,0.975))
summat[8,] <- quantile(sqrt(unlist(sigma2xes)),probs=c(0.025,0.25,0.5,0.75,0.975))
summat[9,] <- summary(corrx)$quantiles
summat[10:12,] <- summary(betaee)$quantiles
summat[13:15,] <- summary(betaes)$quantiles

summat <- data.frame(summat)
names(summat) <- c("2.5\\%","25\\%","50\\%","75\\%","97.5\\%")
row.names(summat) <- c("$\\mu_x^{EE}$","$\\mu_x^{\\Delta ES}$","$\\sigma_{yee}$","$\\sigma_{yes}$",
                       "$\\sigma_{wee}$","$\\sigma_{wes}$","$\\sigma_{xee}$","$\\sigma_{xes}$",
                       "$\\rho$","$\\gamma_{1,ee}$","$\\gamma_{2,ee}$","$\\gamma_{3,ee}$",
                       "$\\gamma_{1,es}$","$\\gamma_{2,es}$","$\\gamma_{3,es}$")

print(xtable(summat), sanitize.text.function=function(x){x})

#---------------------------------------------------------
#plots of spline function
#ee
ind <- sample(1:nrow(chain1$meanfcnee),500)
y <- melt(t(chain1$meanfcnee)[,ind])
x <- melt(t(chain1$latentxee)[,ind])
datee <- cbind(x,y$value)
names(datee) = c("x1","variable", "xval", "yval")

#cs95ee <- data.frame(t(apply(chain1$meanfcnee,2,quantile,probs=c(0.025,0.975))),x=rowMeans(t(chain1$latentxee)))
#names(cs95ee) <- c("lower","upper","x")

p1 <- ggplot()  + geom_line(data = datee,aes(x = xval, y = yval, group = variable),alpha=0.02) +
  geom_line(aes(x=rowMeans(t(chain1$latentxee)),y=rowMeans(t(chain1$meanfcnee))),col="red",size=1) + 
  #geom_line(aes(x=simdata$xee,y=eecurve(simdata$xee)),col="blue",size=1)+
  #geom_ribbon(data=cs95ee,aes(x=x,ymin=lower,ymax=upper,linetype=NA),colour="blue",alpha=0.3) +
  geom_point(aes(x=simdata$xee,y=yeeb)) + xlab("Truth") + ylab("Observed") + theme_bw()

#es
ind <- sample(1:nrow(chain1$meanfcnes),500)
y <- melt(t(chain1$meanfcnes)[,ind])
x <- melt(t(chain1$latentxes)[,ind])
dates <- cbind(x,y$value)
names(dates) = c("x1","variable", "xval", "yval")

#cs95es <- data.frame(t(apply(chain1$meanfcnes,2,quantile,probs=c(0.025,0.975))),x=rowMeans(t(chain1$latentxes)))
#names(cs95es) <- c("lower","upper","x")

p2 <- ggplot() +   geom_line(data = dates,aes(x = xval, y = yval, group = variable),alpha=0.04) +
  geom_line(aes(x=rowMeans(t(chain1$latentxes)),y=rowMeans(t(chain1$meanfcnes))),col="red",size=1) + 
  #geom_line(aes(x=simdata$xee,y=eecurve(simdata$xee)),col="blue",size=1)+
  #geom_ribbon(data=cs95es,aes(x=x,ymin=lower,ymax=upper,linetype=NA),colour="blue",alpha=0.3) +
  geom_point(aes(x=simdata$xes,y=yesb)) + xlab("Truth") + ylab("Observed") + theme_bw()

grid.arrange(p1,p2,nrow=1)

#---------------------------------------------------
#calibration
a <- 100 # 17

ypredee1 <-  yee[a,1] #2500
ypredee2 <-  yee[296,1] #3200
ypredee3 <-  yee[184,1] #4000

ypredes1 <-  yes[a,1]#-250
ypredes2 <-  yes[296,1]#320
ypredes3 <-  yes[184,1]#90

xtrueee1 <- xee[a]
xtrueee2 <- xee[296]
xtrueee3 <- xee[184]

xtruees1 <- xes[a]
xtruees2 <- xes[296]
xtruees3 <- xes[184]

dem1 <- Z[a,] #c(0,18,25)
dem2 <- Z[296,]#c(1,32,20)
dem3 <- Z[184,]#c(1,27,34)
nr <- 1000

ng <- apply(as.matrix(gammaee),1,function(x) sum(x!=0))
nges <- apply(as.matrix(gammaes),1,function(x) sum(x!=0))
nk <- as.numeric(as.matrix(kee))
nkes <- as.numeric(as.matrix(kes))

cee1 <- callibrate(ypredee1,dem1,as.matrix(latentxee),as.matrix(ree),as.matrix(betaee),as.matrix(gammaee),nk,ng,nr,min=1000,max=4500)
#cee2 <- callibrate(ypredee2,dem2,as.matrix(latentxee),as.matrix(ree),as.matrix(betaee),as.matrix(gammaee),nk,ng,nr,min=1500,max=4500)
#cee3 <- callibrate(ypredee3,dem3,as.matrix(latentxee),as.matrix(ree),as.matrix(betaee),as.matrix(gammaee),nk,ng,nr,min=1500,max=4500)

ceei <- data.frame(t(apply(cbind(cee1,cee2,cee3),2,quantile,probs=c(0.025,0.5,0.975))))
names(ceei) <- c("Lower","Median","Upper")
ceei$Observed <- c(ypredee1,ypredee2,ypredee3)

ces1 <- callibrate(ypredes1,dem1,as.matrix(latentxes),as.matrix(res),as.matrix(betaee),as.matrix(gammaes),nkes,nges,nr,min=-300,max=500)
#ces2 <- callibrate(ypredes2,dem2,as.matrix(latentxes),as.matrix(res),as.matrix(betaee),as.matrix(gammaes),nkes,nges,nr,min=-100,max=100)
#ces3 <- callibrate(ypredes3,dem3,as.matrix(latentxes),as.matrix(res),as.matrix(betaes),as.matrix(gammaes),nkes,nges,nr,min=-500,max=300)

cesi <- data.frame(t(apply(cbind(ces1,ces2,ces3),2,quantile,probs=c(0.025,0.5,0.975))))
names(cesi) <- c("Lower","Median","Upper")
cesi$Observed <- c(ypredes1,ypredes2,ypredes3)

# c1 <- ggplot() +   geom_line(data = datee,aes(x = (xval), y = yval, group = variable),alpha=0.04) +
#      geom_segment(data=ceei,aes(x=Lower,xend=Upper,y=Observed,yend=Observed),colour="blue",size=1) +
#      geom_point(data=ceei,aes(x=Median,y=Observed),colour="red",size=2) + xlab("Truth") + ylab("Observed") + theme_bw()
# 
# c2 <- ggplot() +   geom_line(data = dates,aes(x = (xval), y = yval, group = variable),alpha=0.04) +
#   geom_segment(data=cesi,aes(x=Lower,xend=Upper,y=Observed,yend=Observed),colour="blue",size=1) +
#   geom_point(data=cesi,aes(x=Median,y=Observed),colour="red",size=2) + xlab("Truth") + ylab("Observed") + theme_bw()
#   #geom_boxplot(aes(x=(cee1),y=(ypredee1))) + geom_boxplot(aes(x=cee2,y=ypredee2)) 
# 
# grid.arrange(c1,c2,nrow=1)

dfcal <- data.frame(cbind(cee1,cee2,cee3,ces1,ces2,ces3))
names(dfcal) <- c("Calibrated EE 1","Calibrated EE 2","Calibrated EE 3","Calibrated ES 1","Calibrated ES 2","Calibrated ES 3")
mdfcal <- melt(dfcal)
mdfcal$obs <- c(rep(ypredee1,nr),rep(ypredee2,nr),rep(ypredee3,nr),rep(ypredes1,nr),rep(ypredes2,nr),rep(ypredes3,nr))
mdfcal$true <- c(rep(xtrueee1,nr),rep(xtrueee2,nr),rep(xtrueee3,nr),rep(xtruees1,nr),rep(xtruees2,nr),rep(xtruees3,nr))

ggplot(data=mdfcal) + geom_histogram(aes(x=value),bins=20) + geom_vline(aes(xintercept=obs),colour="red") + geom_vline(aes(xintercept=true),colour="blue",linetype=2) + facet_wrap(~variable,scales="free") + theme_bw()



names(ceei)[1:3] <- c("2.5%","50%","97.5%")
names(cesi)[1:3] <- c("2.5%","50%","97.5%")
ceei$Truth <- c(xee[100],xee[296],xee[184])
cesi$Truth <- c(xes[100],xes[296],xes[184])
print(xtable(ceei,align="cccc|cc",caption="95\\% credible interval for calibration estimate for cheap EE measurements for Skewed Errors",label="calibratedee"),include.rownames=FALSE)
print(xtable(cesi,align="cccc|cc",caption="95\\% credible interval for calibration estimate for cheap $\\Delta$ES measurements for Skewed Errors",label="calibratedes"),include.rownames=FALSE)


#--------------------------------------------------
#same as above except do it for every observation in yeeb and yesb and give
#observation its predicted interval
nsim <- 10
cee <- matrix(0,ncol=length(yeeb),nrow=nsim)
for(i in 1:length(yeeb)){
  cee[,i] <- callibrate(yeeb[i],t(Z[i,]),as.matrix(latentxee),as.matrix(ree),as.matrix(betaee),as.matrix(gammaee),nk,ng,nsim,min=yeeb[i]-500,max=yeeb[i]+500)
}

ceei <- data.frame(t(apply(cee,2,quantile,probs=c(0.025,0.5,0.975))))
names(ceei) <- c("Lower","Median","Upper")
ceei$Observed <- yeeb
ceei$true <- xee
# 
# ces1 <- callibrate(ypredes1,0,as.matrix(latentxes),as.matrix(res),0,as.matrix(gammaes),nkes,nges,500,min=-200,max=100)
# ces2 <- callibrate(ypredes2,0,as.matrix(latentxes),as.matrix(res),0,as.matrix(gammaes),nkes,nges,500,min=0,max=200)
# ces3 <- callibrate(ypredes3,c(1,27,27),as.matrix(latentxes),as.matrix(res),as.matrix(betaes),as.matrix(gammaes),nkes,nges,500,min=-200,max=300)
# 
# cesi <- data.frame(t(apply(cbind(ces1,ces2,ces3),2,quantile,probs=c(0.025,0.5,0.975))))
# names(cesi) <- c("Lower","Median","Upper")
# cesi$Observed <- c(ypredes1,ypredes2,ypredes3)
# 
c1 <- ggplot() + geom_point(aes(x=xee,y=yeeb)) +
  geom_line(data = datee,aes(x = (xval), y = yval, group = variable),alpha=0.04) +
  geom_segment(data=ceei,aes(x=Lower,xend=Upper,y=Observed,yend=Observed),colour="blue",size=2) +
  geom_point(data=ceei,aes(x=Median,y=Observed),colour="red",size=3) + xlab("Truth") + ylab("Observed") + theme_bw()
# 
# c2 <- ggplot() +   geom_line(data = dates,aes(x = (xval), y = yval, group = variable),alpha=0.04) +
#   geom_segment(data=cesi,aes(x=Lower,xend=Upper,y=Observed,yend=Observed),colour="blue",size=2) +
#   geom_point(data=cesi,aes(x=Median,y=Observed),colour="red",size=3) + xlab("Truth") + ylab("Observed") + theme_bw()
# #geom_boxplot(aes(x=(cee1),y=(ypredee1))) + geom_boxplot(aes(x=cee2,y=ypredee2)) 
# 
# grid.arrange(c1,c2,nrow=1)
# 
#----------------------------------------------#
#extract knot locations

#********still need to incorporate linear terms
callibrate <- function(y,z,latentx,knotlocations,betacoef,gammacoef,nk,ng,ndraws,min=0,max=5500){
  #given a value of y, we want to calculate the corresponding value of x
  # y is value to predict with
  # z is vector of observed error free covariates corresponding to individual of interest
  # latent x is a matrix of latent x variables for all individuals, each col is individual and each 
  # row represents another draw
  # knotlocations is matrix of knot locations
  indicies <- sample(1:nrow(gammacoef),ndraws)
  print(indicies)
  prediction <- rep(0,ndraws)
  if(betacoef==0){
    betacoef <- matrix(0,nrow=nrow(gammacoef),ncol=1)
  }
  ct <- 1
  for(i in indicies){
    basis_matrix <- bs(latentx[i,],knots=knotlocations[i,1:nk[i]],intercept=TRUE)
    coefs <- gammacoef[i,1:ng[i]]
    my_inverse <- inverse(function(x) {return(coefs%*%(as.matrix(predict(basis_matrix,x))[1,]))},min,max)
    prediction[ct] <- my_inverse(y) + t(z)%*%betacoef[i,]
    ct <- ct+1
  }
  return(prediction)
}

ng <- apply(as.matrix(gammaee),1,function(x) sum(x!=0))
nges <- apply(as.matrix(gammaes),1,function(x) sum(x!=0))
nk <- as.numeric(as.matrix(kee))
nkes <- as.numeric(as.matrix(kes))


a1 <- callibrate(2900,0,as.matrix(latentxee),as.matrix(ree),0,as.matrix(gammaee),nk,ng,1000,min=0,max=4500)
hist(a1)
a2 <- callibrate(500,0,as.matrix(latentxes),as.matrix(res),0,as.matrix(gammaes),nkes,nges,300,min=-400,max=400)

