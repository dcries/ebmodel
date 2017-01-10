nrep <- 500
yeesd <- rep(0,nrep)
yessd <- rep(0,nrep)

for(k in 1:nrep){
n <- 300
#known covariates
zg <- rbinom(n,1,0.5) #gender indicator
zb <- rnorm(n,27,5) #bmi
za <- runif(n,20,40) #age
#ze <- rbinom(n,1,)  #education
#zr <-  #race

geg <- params[3]#200 #indicator level for male for ee
geb <- params[4]#3.5 #slope for bmi for ee
gea <- params[5]#-2.5 #slope for age for ee

gig <- params[6]#-200 #indicator level for male for ei
gib <-params[7]# 4 #slope for bmi for ei
gia <- params[8]#-3 #slope for age for ei


#slopes for true ee/ei
#note zig and gig are much different, zig refers to ei and gig to es, a mistake
zeg <- 400 #indicator level for male for ee
zeb <- 10 #slope for bmi for ee
zea <- -7 #slope for age for ee

zig <- 350 #indicator level for male for ei
zib <- 14 #slope for bmi for ei
zia <- -6 #slope for age for ei

#mean daily EI
#mean daily EI
mei1 <- 1500 + zig*zg + zib*zb + zia*za
mei2 <- 1900 + zig*zg + zib*zb + zia*za
mei3 <- 2100 + zig*zg + zib*zb + zia*za
mei4 <- 2900 + zig*zg + zib*zb + zia*za
mei5 <- 3200 + zig*zg + zib*zb + zia*za
#sd EI
sdei1 <- 50 +0#50
sdei2 <- 80 +0#50
sdei3 <- 60 +0#90
sdei4 <- 400 +0#40
sdei5 <- 220 +0#90

#mean daily EE = mean daily EI
mee1 <- 1500 + 130 + zeg*zg + zeb*zb + zea*za # 1500
mee2 <- 1900 + 130 + zeg*zg + zeb*zb + zea*za # 2000
mee3 <- 2100 + 80 + zeg*zg + zeb*zb + zea*za # 2300
mee4 <- 2900 + 130 + zeg*zg + zeb*zb + zea*za # 2600
mee5 <- 3200 + 130 + zeg*zg + zeb*zb + zea*za #- 3500
mee <- c(mee1,mee2,mee3,mee4,mee5)
#sd EE
sdee1 <- sdei1 + 20
sdee2 <- sdei2 + 15
sdee3 <- sdei3 + 30
sdee4 <- sdei4 + 10
sdee5 <- sdei5 + 30
#cor btwn the two
rho <- 0.4376205 #cor(eb$Energy_expenditure,eb$energy,use="complete.obs")
#cor btwn ee and es, calculate
rhos <- -(1-rho)
#sd es, calculated
#sdes <- sqrt(sdee^2+sdei^2-2*sdee*sdei*rho)
#cov matrix btwn the two
xcov1 <- matrix(c(sdei1^2, rho*sdei1*sdee1 , rho*sdei1*sdee1, sdee1^2),ncol=2,byrow=T)
xcov2 <- matrix(c(sdei2^2, rho*sdei2*sdee2 , rho*sdei2*sdee2, sdee2^2),ncol=2,byrow=T)
xcov3 <- matrix(c(sdei3^2, rho*sdei3*sdee3 , rho*sdei3*sdee3, sdee3^2),ncol=2,byrow=T)
xcov4 <- matrix(c(sdei4^2, rho*sdei4*sdee4 , rho*sdei4*sdee4, sdee4^2),ncol=2,byrow=T)
xcov5 <- matrix(c(sdei5^2, rho*sdei5*sdee5 , rho*sdei5*sdee5, sdee5^2),ncol=2,byrow=T)

#xcovs <- matrix(c(sdee^2, rhos*sdes*sdee , rhos*sdes*sdee, sdes^2),ncol=2,byrow=T)
#intercept of yee bias
be0 <- params[1]#100
#intercept of yes bias
bs0 <- params[2]#50
#slope of yee bias
#be1 <- 0.6
#slopt of yes bias
#bs1 <- 1.2
#simulate n observations of daily EI and EE, mixture
x1 <- matrix(c(rep(0,n/5),rep(0,n/5)),ncol=2)
x2 <- matrix(c(rep(0,n/5),rep(0,n/5)),ncol=2)
x3 <- matrix(c(rep(0,n/5),rep(0,n/5)),ncol=2)
x4 <- matrix(c(rep(0,n/5),rep(0,n/5)),ncol=2)
x5 <- matrix(c(rep(0,n/5),rep(0,n/5)),ncol=2)

# for(i in 1:(n/5)){
#   x1[i,] <- mvrnorm(1,c(mei1[i],mee1[i]),xcov1)
#   x2[i,] <- mvrnorm(1,c(mei2[(i+n/5)],mee2[(i+n/5)]),xcov2)
#   x3[i,] <- mvrnorm(1,c(mei3[(i+2*n/5)],mee3[(i+2*n/5)]),xcov3)
#   x4[i,] <- mvrnorm(1,c(mei4[(i+3*n/5)],mee4[(i+3*n/5)]),xcov4)
#   x5[i,] <- mvrnorm(1,c(mei5[(i+4*n/5)],mee5[(i+4*n/5)]),xcov5)
# }
# 
df <- 4
for(i in 1:(n/5)){
  x1[i,] <- rmvt(1,delta=c(mei1[i],mee1[i]),sigma=xcov1*(df-2)/df,df=df)
  x2[i,] <- rmvt(1,delta=c(mei2[(i+n/5)],mee2[(i+n/5)]),sigma=xcov2*(df-2)/df,df=df)
  x3[i,] <- rmvt(1,delta=c(mei3[(i+2*n/5)],mee3[(i+2*n/5)]),sigma=xcov3*(df-2)/df,df=df)
  x4[i,] <- rmvt(1,delta=c(mei4[(i+3*n/5)],mee4[(i+3*n/5)]),sigma=xcov4*(df-2)/df,df=df)
  x5[i,] <- rmvt(1,delta=c(mei5[(i+4*n/5)],mee5[(i+4*n/5)]),sigma=xcov5*(df-2)/df,df=df)
}

# x1 <- mvrnorm(n/5,c(mei1,mee1),xcov1)
# x2 <- mvrnorm(n/5,c(mei2,mee2),xcov2)
# x3 <- mvrnorm(n/5,c(mei3,mee3),xcov3)
# x4 <- mvrnorm(n/5,c(mei4,mee4),xcov4)
# x5 <- mvrnorm(n/5,c(mei5,mee5),xcov5)

#for bimodal error model
mix <- rbinom(n,1,0.5) 

xei <- c(x1[,1],x2[,1],x3[,1],x4[,1],x5[,1])
xee <- c(x1[,2],x2[,2],x3[,2],x4[,2],x5[,2])

#within person variability
dmatrix <- matrix(c(50^2,-2000,-2000,150^2),ncol=2,byrow=TRUE) #EI first row
#delta <- list()
#for(j in 1:nrep){
  delta = mvrnorm(n,rep(0,2),dmatrix)
#}


#calculate delta ES, positive change => ei > ee
#xes <- xei - xee

xes <- sample(c(rnorm(n/5,-38,23),rnorm(n/5,38,23),rnorm(n/5,0,34),rnorm(n/5,-80,85),rnorm(n/5,80,85)),n,replace=FALSE)

#calculate true values T_ij of EE,EI, ES
#   tei1 <- xei + delta1[,1]
#   tei2 <- xei + delta2[,1]
#   tee1 <- xee + delta1[,2]
#   tee2 <- xee + delta2[,2]
#   tes1 <- tei1 - tee1
#   tes2 <- tei2 - tee2
#calculate DLW est of EE



#yee <- be0+eecurve(xee+ delta[,2],4000,2100,0.002)+geg*zg+geb*zb+gea*za + mix*rnorm(n,-350,sqrt(380^2-350^2)) + (1-mix)*rnorm(n,350,sqrt(380^2-350^2)) 
#calc cheap ES
#yes <- bs0+escurve(xes+ delta[,1],800,k=0.07)+gig*zg+gib*zb+gia*za + mix*rnorm(n,-190,sqrt(210^2-190^2)) + (1-mix)*rnorm(n,190,sqrt(210^2-190^2)) 
yee <- be0+eecurve(xee+ delta[,2])+geg*zg+geb*zb+gea*za + rsnorm(n,0,380,10) 
#calc cheap ES
yes <- bs0+escurve(xes+ delta[,1],k=0.1)+gig*zg+gib*zb+gia*za + rsnorm(n,0,210,10) 

     # yee <- be0+eecurve(xee+ delta[,2],4000,2100,0.002)+geg*zg+geb*zb+gea*za + rnorm(n,0,380) #truth 408.534
    #  yes <- bs0+escurve(xes+ delta[,1],800,k=0.07)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) #truth 216.1748
      yeesd[k] <- sd(yee-(be0+eecurve(xee,4000,2100,0.002)+geg*zg+geb*zb+gea*za))
      yessd[k] <- sd(yes-(bs0+escurve(xes,k=0.1)+gig*zg+gib*zb+gia*za))
      #print(k)
  }