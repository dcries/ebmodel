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
Rcpp::sourceCpp('C:\\Users\\dcries\\github\\ebmodel\\bivar_bvnmcmc1_qp.cpp')

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
burn <- 500
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
currentxee <- wee[,1];currentxee2 <- yee[,2];currentxee3 <- wee[,2]+rnorm(n,0,300)
currentxes <- wes[,1];currentxes2 <- yes[,1];currentxes3 <- wes[,2]+rnorm(n,0,100)
currentmuee <- rep(2600,h);currentmuee2 <- rep(1800,h);currentmuee3 <- rep(3600,h)#runif(h,1800,3800)
currentmues <- rep(0,h);currentmues2 <- rep(100,h);currentmues3 <- rep(-100,h)#runif(h,-400,400)
#currentmuee <- 2600#runif(h,1800,3800)
#currentmues <- 0#runif(h,-400,400)
currentpi <- rep(1/h,h)
currentzeta <- sample(1:h,n,replace=T)
currentv <- rep(0.3,h)
currentpredee <- wee[,1]
currentpredes <- wes[,1]
currentsigma2ee <- 400^2;currentsigma2ee2 <- 200^2;currentsigma2ee3 <- 600^2
currentsigma2ve <- 250^2;currentsigma2ve2 <- 150^2;currentsigma2ve3 <- 350^2
currentsigma2es <- 220^2;currentsigma2es2 <- 120^2;currentsigma2es3 <- 320^2
currentsigma2vs <- 80^2;currentsigma2vs2 <- 180^2;currentsigma2vs3 <- 30^2
currentalpha <- 1;currentalpha2 <- 4;currentalpha3 <- 0.3
tunevar = c(100^2,25^2)#c(100^2,25^2)
tunecor = -0.2
currentsigma2x <- c(400^2,400^2);currentsigma2x2 <- c(800^2,800^2);currentsigma2x3 <- c(100^2,100^2)
currentbetaee <- matrix(c(0,0,0),nrow=1);currentbetaee2 <- matrix(c(-50,-15,-15),nrow=1);currentbetaee3 <- matrix(c(100,15,15),nrow=1)
currentbetaes <- matrix(c(0,0,0),nrow=1);currentbetaes2 <- matrix(c(-100,-15,-15),nrow=1);currentbetaes3 <- matrix(c(100,15,15),nrow=1)


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
initial <- list(currentkee=currentkee,currentkes=currentkes,ck=ck,knotsee=knotsee,
                knotses=knotses,currentxee=currentxee,currentxes=currentxes,
                currentv=currentv,currentpi=currentpi,currentalpha=currentalpha,
                currentzeta=currentzeta,currentpredee=currentpredee,currentpredes=currentpredes,
                currentmuee=currentmuee,currentmues=currentmues,currentsigma2ee=currentsigma2ee,
                currentsigma2es=currentsigma2es,currentsigma2ve=currentsigma2ve,
                currentsigma2vs=currentsigma2vs,currentsigma2x=currentsigma2x,
                tunevar=tunevar,currentbetaee=currentbetaee,currentbetaes=currentbetaes,
                tunecor=tunecor)


initial2 <- list(currentkee=currentkee2,currentkes=currentkes2,ck=ck,knotsee=knotsee2,
                 knotses=knotses2,currentxee=currentxee2,currentxes=currentxes2,
                 currentv=currentv,currentpi=currentpi,currentalpha=currentalpha2,
                 currentzeta=currentzeta,currentpredee=currentpredee,currentpredes=currentpredes,
                 currentmuee=currentmuee2,currentmues=currentmues2,currentsigma2ee=currentsigma2ee2,
                 currentsigma2es=currentsigma2es2,currentsigma2ve=currentsigma2ve2,
                 currentsigma2vs=currentsigma2vs2,currentsigma2x=currentsigma2x2,
                 tunevar=tunevar,currentbetaee=currentbetaee2,currentbetaes=currentbetaes2,
                 tunecor=tunecor)


initial3 <- list(currentkee=currentkee3,currentkes=currentkes3,ck=ck,knotsee=knotsee3,
                 knotses=knotses3,currentxee=currentxee3,currentxes=currentxes3,
                 currentv=currentv,currentpi=currentpi,currentalpha=currentalpha3,
                 currentzeta=currentzeta,currentpredee=currentpredee,currentpredes=currentpredes,
                 currentmuee=currentmuee3,currentmues=currentmues3,currentsigma2ee=currentsigma2ee3,
                 currentsigma2es=currentsigma2es3,currentsigma2ve=currentsigma2ve3,
                 currentsigma2vs=currentsigma2vs3,currentsigma2x=currentsigma2x3,
                 tunevar=tunevar,currentbetaee=currentbetaee3,currentbetaes=currentbetaes3,
                 tunecor=tunecor)

prior <- list(lambda=lambda,ae=ae,be=be,av=av,bv=bv,a_alp=a_alp,
              b_alp=b_alp,d=d,m=m,v2=v2,psi=psi,Vb=Vb,Mb=Mb)

#chain1=mcmc_bvn(yee,yes,wee,wes,Z,initial,prior,nreps,burn,maxkt,my_bs)
#chain2=mcmc_bvn(yee,yes,wee,wes,Z,initial2,prior,nreps,burn,maxkt,my_bs)
#chain3=mcmc_bvn(yee,yes,wee,wes,Z,initial3,prior,nreps,burn,maxkt,my_bs)

chain1=mcmc_bvn_qp(yee,yes,wee,wes,Z,initial,prior,nreps,burn,maxkt,my_bs,my_qp)
chain2=mcmc_bvn_qp(yee,yes,wee,wes,Z,initial2,prior,nreps,burn,maxkt,my_bs,my_qp)
chain3=mcmc_bvn_qp(yee,yes,wee,wes,Z,initial3,prior,nreps,burn,maxkt,my_bs,my_qp)

#---------------------------------------------------------------------------------#
#combine lists
#mapply(c, a1, a2, SIMPLIFY=FALSE)

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
kee <- mcmc.list(mcmc(chain1$kee),mcmc(chain2$kee),mcmc(chain3$kee))
kes <- mcmc.list(mcmc(chain1$kes),mcmc(chain2$kes),mcmc(chain3$kes))
betaee <- mcmc.list(mcmc(chain1$betaee),mcmc(chain2$betaee),mcmc(chain3$betaee))
betaes <- mcmc.list(mcmc(chain1$betaes),mcmc(chain2$betaes),mcmc(chain3$betaes))

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

out3x <- list(chain1=chain1,chain2=chain2,chain3=chain3,ppan=ppan)
save(out3x,file="//my.files.iastate.edu/Users/dcries/Desktop/bvn_mcmc3x.RData")
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


#ee
ind <- sample(1:nrow(chain1$meanfcnee),100)
y <- melt(t(chain1$meanfcnee)[,ind])
x <- melt(t(chain1$latentxee)[,ind])
dat <- cbind(x,y$value)
names(dat) = c("x1","variable", "xval", "yval")

#cs95ee <- data.frame(t(apply(chain1$meanfcnee,2,quantile,probs=c(0.025,0.975))),x=rowMeans(t(chain1$latentxee)))
#names(cs95ee) <- c("lower","upper","x")

p1 <- ggplot()  + geom_line(data = dat,aes(x = xval, y = yval, group = variable),alpha=0.02) +
  geom_line(aes(x=rowMeans(t(chain1$latentxee)),y=rowMeans(t(chain1$meanfcnee))),col="red",size=1) + 
  #geom_line(aes(x=simdata$xee,y=eecurve(simdata$xee)),col="blue",size=1)+
  #geom_ribbon(data=cs95ee,aes(x=x,ymin=lower,ymax=upper,linetype=NA),colour="blue",alpha=0.3) +
  geom_point(aes(x=simdata$xee,y=yeeb)) + xlab("Truth") + ylab("Observed") + theme_bw()

#es
ind <- sample(1:nrow(chain1$meanfcnes),100)
y <- melt(t(chain1$meanfcnes)[,ind])
x <- melt(t(chain1$latentxes)[,ind])
dat <- cbind(x,y$value)
names(dat) = c("x1","variable", "xval", "yval")

#cs95es <- data.frame(t(apply(chain1$meanfcnes,2,quantile,probs=c(0.025,0.975))),x=rowMeans(t(chain1$latentxes)))
#names(cs95es) <- c("lower","upper","x")

p2 <- ggplot() +   geom_line(data = dat,aes(x = xval, y = yval, group = variable),alpha=0.04) +
  geom_line(aes(x=rowMeans(t(chain1$latentxes)),y=rowMeans(t(chain1$meanfcnes))),col="red",size=1) + 
  #geom_line(aes(x=simdata$xee,y=eecurve(simdata$xee)),col="blue",size=1)+
  #geom_ribbon(data=cs95es,aes(x=x,ymin=lower,ymax=upper,linetype=NA),colour="blue",alpha=0.3) +
  geom_point(aes(x=simdata$xes,y=yesb)) + xlab("Truth") + ylab("Observed") + theme_bw()

grid.arrange(p1,p2,nrow=1)



#----------------------------------------------#
#extract knot locations
ree <- chain1$ree
gammaee <- chain1$gammaee
apply(ree,1,function(x) length(unique(x)))-1 #gives kee
apply(gammaee,1,function(x) sum(x!=0)) #gives kee

