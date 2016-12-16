
d <- matrix(0,ncol=3,nrow=500)
for(i in 1:500){
  simdata <- generate_data5(params,dist=1,nrep=2)

  yeeb <- rowMeans(simdata$yee)
  xee <- rowMeans(simdata$wee)
  zg <- simdata$zg #<- rbinom(n,1,0.5) #gender indicator
  zb <- simdata$zb#<- rnorm(n,27,5) #bmi
  za <- simdata$za#<- runif(n,20,40) #age
  
  Z= cbind(zg,zb,za)
  

  m1 <- lm(yeeb~Z+bs(xee,df=6,intercept=F))
  d[i,] <- coef(m1)[2:4]
}

summary(d)


d <- matrix(0,ncol=3,nrow=500)
for(i in 1:500){
  simdata <- generate_data5(params,dist=1,nrep=2)
  
  yeeb <- rowMeans(simdata$yee)
  xee <- rowMeans(simdata$wee)
  zg <- simdata$zg #<- rbinom(n,1,0.5) #gender indicator
  zb <- simdata$zb#<- rnorm(n,27,5) #bmi
  za <- simdata$za#<- runif(n,20,40) #age
  
  Z= cbind(zg,zb,za)
  
  
  m1 <- lm(yeeb~bs(xee,df=6,intercept=T))
  m2 <- lm(yeeb-predict(m1)~Z+0)
 
  d[i,] <- coef(m2)[1:3]
}

summary(d)