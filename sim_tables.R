library(xtable)
setwd("C:/Users/dcries/github/ebmodel/workspace")
load("bvn_2_1.RData")
load("bvn_2_2.RData")
load("bvn_2_3.RData")
load("bvn_4_1.RData")
load("bvn_4_2.RData")
load("bvn_4_3.RData")

bvn_4_2 <- bvn_2_4
bvn_2_3 <- bvn_3_2
bvn_4_3 <- bvn_3_4

load("dpmm_4_1.RData")
load("dpmm_4_2.RData")
load("dpmm_4_3.RData")

#note each of these loads to object named out
load("simple_2_1.RData")
simple_2_1 <- out
load("simple_2_2.RData")
simple_2_2 <- out
load("simple_2_3.RData")
simple_2_3 <- out
load("simple_4_1.RData")
simple_4_1 <- out
load("simple_4_2.RData")
simple_4_2 <- out
load("simple_4_3.RData")
simple_4_3 <- out


#-----------------------------------------------
#for bvn only
out <- bvn_2_1
pmse <- out$pmse
model <- out$model[,-c(1,2,7,8,9)]
indcheck <- out$indcheck

#----------------------------------------------
#for dpmm only
out <- dpmm_4_1
pmse <- out$pmse
model <- out$model
indcheck <- out$indcheck

#----------------------------------------------
#for simple only



#order sigmaeee,sigmaees,sigmavee,sigmaves,betaee,betaes
true <- c(0,0,250,72.862,300,14,-7,-200,8,-5)


meanest <- colMeans(model)
rmse <- apply(model,2,sd)
bias <- meanest-true

df <- data.frame(True=true,MeanEst=meanest,StdErr=rmse,Bias=bias)
df <- t(df)
df <- as.data.frame(df)
names(df) <- c("$\\sigma_{yee}$","$\\sigma_{yes}$","$\\sigma_{wee}$","$\\sigma_{wes}$",
               "$\\gamma_{1,ee}$","$\\gamma_{2,ee}$","$\\gamma_{3,ee}$","$\\gamma_{1,es}$",
               "$\\gamma_{2,es}$","$\\gamma_{3,es}$")

print(xtable(df),sanitize.text.function = function(x){x})

#---------------------------------------
#pmse ee matrix
pmatee <- matrix(0,nrow=6,ncol=4)

pmatee[1,1] <- mean(simple_2_1$pmse0[,1])
pmatee[2,1] <- mean(simple_4_1$pmse0[,1])
pmatee[3,1] <- mean(simple_2_2$pmse0[,1])
pmatee[4,1] <- mean(simple_4_2$pmse0[,1])
pmatee[5,1] <- mean(simple_2_3$pmse0[,1])
pmatee[6,1] <- mean(simple_4_3$pmse0[,1])

pmatee[1,2] <- mean(simple_2_1$pmse1[,1])
pmatee[2,2] <- mean(simple_4_1$pmse1[,1])
pmatee[3,2] <- mean(simple_2_2$pmse1[,1])
pmatee[4,2] <- mean(simple_4_2$pmse1[,1])
pmatee[5,2] <- mean(simple_2_3$pmse1[,1])
pmatee[6,2] <- mean(simple_4_3$pmse1[,1])

pmatee[1,3] <- mean(bvn_2_1$pmse[,1])
pmatee[2,3] <- mean(bvn_4_1$pmse[,1])
pmatee[3,3] <- mean(bvn_2_2$pmse[,1])
pmatee[4,3] <- mean(bvn_4_2$pmse[,1])
pmatee[5,3] <- mean(bvn_2_3$pmse[,1])
pmatee[6,3] <- mean(bvn_4_3$pmse[,1])

pmatee[1,4] <- 0
pmatee[2,4] <- mean(dpmm_4_1$pmse[,1])
pmatee[3,4] <- 0
pmatee[4,4] <- mean(dpmm_4_2$pmse[,1])
pmatee[5,4] <- 0
pmatee[6,4] <- mean(dpmm_4_3$pmse[,1])
#---------------------------------------

pmates <- matrix(0,nrow=6,ncol=4)
pmates[1,1] <- mean(simple_2_1$pmse0[,2])
pmates[2,1] <- mean(simple_4_1$pmse0[,2])
pmates[3,1] <- mean(simple_2_2$pmse0[,2])
pmates[4,1] <- mean(simple_4_2$pmse0[,2])
pmates[5,1] <- mean(simple_2_3$pmse0[,2])
pmates[6,1] <- mean(simple_4_3$pmse0[,2])

pmates[1,2] <- mean(simple_2_1$pmse1[,2])
pmates[2,2] <- mean(simple_4_1$pmse1[,2])
pmates[3,2] <- mean(simple_2_2$pmse1[,2])
pmates[4,2] <- mean(simple_4_2$pmse1[,2])
pmates[5,2] <- mean(simple_2_3$pmse1[,2])
pmates[6,2] <- mean(simple_4_3$pmse1[,2])

pmates[1,3] <- mean(bvn_2_1$pmse[,2])
pmates[2,3] <- mean(bvn_4_1$pmse[,2])
pmates[3,3] <- mean(bvn_2_2$pmse[,2])
pmates[4,3] <- mean(bvn_4_2$pmse[,2])
pmates[5,3] <- mean(bvn_2_3$pmse[,2])
pmates[6,3] <- mean(bvn_4_3$pmse[,2])

pmates[1,4] <- 0
pmates[2,4] <- mean(dpmm_4_1$pmse[,2])
pmates[3,4] <- 0
pmates[4,4] <- mean(dpmm_4_2$pmse[,2])
pmates[5,4] <- 0
pmates[6,4] <- mean(dpmm_4_3$pmse[,2])

#----------------------------------------------#

#true20 <- c(420,420,344,344,-7,-7,14,14,300,300,-5,-5,8,8,-200,-200)
true20 <- c(405.5,405.5,344,344,300,300,14,14,-7,-7,-200,-200,8,8,-5,-5)

true2 <- c(405.5,405.5,344,344,250,250,72.862,72.862,300,300,14,14,-7,-7,-200,-200,8,8,-5,-5)
true2d <- c(405.5,344,250,72.862,300,14,-7,-200,8,-5)


#----------------------------------------#
#same but for normal data

norm0 <- matrix(0,nrow=16,ncol=4)
norm1 <- matrix(0,nrow=20,ncol=4)
normbvn <- matrix(0,nrow=20,ncol=4)
normdpmm <- matrix(0,nrow=10,ncol=4)


norm0[1,1] <- mean(simple_2_1$model0[,11]);norm0[1,2] <- sd(simple_2_1$model0[,11])
norm0[2,1] <- mean(simple_4_1$model0[,11]);norm0[2,2] <- sd(simple_4_1$model0[,11])
norm0[3,1] <- mean(simple_2_1$model0[,12]);norm0[3,2] <- sd(simple_2_1$model0[,12])
norm0[4,1] <- mean(simple_4_1$model0[,12]);norm0[4,2] <- sd(simple_4_1$model0[,12])
norm0[9,1] <- mean(simple_2_1$model0[,5]);norm0[9,2] <- sd(simple_2_1$model0[,5])
norm0[10,1] <- mean(simple_4_1$model0[,5]);norm0[10,2] <- sd(simple_4_1$model0[,5])
norm0[7,1] <- mean(simple_2_1$model0[,6]);norm0[7,2] <- sd(simple_2_1$model0[,6])
norm0[8,1] <- mean(simple_4_1$model0[,6]);norm0[8,2] <- sd(simple_4_1$model0[,6])
norm0[5,1] <- mean(simple_2_1$model0[,7]);norm0[5,2] <- sd(simple_2_1$model0[,7])
norm0[6,1] <- mean(simple_4_1$model0[,7]);norm0[6,2] <- sd(simple_4_1$model0[,7])
norm0[15,1] <- mean(simple_2_1$model0[,8]);norm0[15,2] <- sd(simple_2_1$model0[,8])  
norm0[16,1] <- mean(simple_4_1$model0[,8]);norm0[16,2] <- sd(simple_4_1$model0[,8])
norm0[13,1] <- mean(simple_2_1$model0[,9]);norm0[13,2] <- sd(simple_2_1$model0[,9])
norm0[14,1] <- mean(simple_4_1$model0[,9]);norm0[14,2] <- sd(simple_4_1$model0[,9])
norm0[11,1] <- mean(simple_2_1$model0[,10]);norm0[11,2] <- sd(simple_2_1$model0[,10])
norm0[12,1] <- mean(simple_4_1$model0[,10]);norm0[12,2] <- sd(simple_4_1$model0[,10])


norm1[1,1] <- mean(simple_2_1$model1[,7]);norm1[1,2] <- sd(simple_2_1$model1[,7])
norm1[2,1] <- mean(simple_4_1$model1[,7]);norm1[2,2] <- sd(simple_4_1$model1[,7])
norm1[3,1] <- mean(simple_2_1$model1[,8]);norm1[3,2] <- sd(simple_2_1$model1[,8])
norm1[4,1] <- mean(simple_4_1$model1[,8]);norm1[4,2] <- sd(simple_4_1$model1[,8])
norm1[5,1] <- mean(simple_2_1$model1[,9]);norm1[5,2] <- sd(simple_2_1$model1[,9])
norm1[6,1] <- mean(simple_4_1$model1[,9]);norm1[6,2] <- sd(simple_4_1$model1[,9])
norm1[7,1] <- mean(simple_2_1$model1[,10]);norm1[7,2] <- sd(simple_2_1$model1[,10])
norm1[8,1] <- mean(simple_4_1$model1[,10]);norm1[8,2] <- sd(simple_4_1$model1[,10])
norm1[9,1] <- mean(simple_2_1$model1[,14]);norm1[9,2] <- sd(simple_2_1$model1[,14])
norm1[10,1] <- mean(simple_4_1$model1[,14]);norm1[10,2] <- sd(simple_4_1$model1[,14])
norm1[11,1] <- mean(simple_2_1$model1[,15]);norm1[11,2] <- sd(simple_2_1$model1[,15])
norm1[12,1] <- mean(simple_4_1$model1[,15]);norm1[12,2] <- sd(simple_4_1$model1[,15])
norm1[13,1] <- mean(simple_2_1$model1[,16]);norm1[13,2] <- sd(simple_2_1$model1[,16])
norm1[14,1] <- mean(simple_4_1$model1[,16]);norm1[14,2] <- sd(simple_4_1$model1[,16])
norm1[15,1] <- mean(simple_2_1$model1[,17]);norm1[15,2] <- sd(simple_2_1$model1[,17])  
norm1[16,1] <- mean(simple_4_1$model1[,17]);norm1[16,2] <- sd(simple_4_1$model1[,17])
norm1[17,1] <- mean(simple_2_1$model1[,18]);norm1[17,2] <- sd(simple_2_1$model1[,18])
norm1[18,1] <- mean(simple_4_1$model1[,18]);norm1[18,2] <- sd(simple_4_1$model1[,18])
norm1[19,1] <- mean(simple_2_1$model1[,19]);norm1[19,2] <- sd(simple_2_1$model1[,19])
norm1[20,1] <- mean(simple_4_1$model1[,19]);norm1[20,2] <- sd(simple_4_1$model1[,19])

normbvn[1,1] <- mean(bvn_2_1$model[,3]);normbvn[1,2] <- sd(bvn_2_1$model[,3])
normbvn[2,1] <- mean(bvn_4_1$model[,3]);normbvn[2,2] <- sd(bvn_4_1$model[,3])
normbvn[3,1] <- mean(bvn_2_1$model[,4]);normbvn[3,2] <- sd(bvn_2_1$model[,4])
normbvn[4,1] <- mean(bvn_4_1$model[,4]);normbvn[4,2] <- sd(bvn_4_1$model[,4])
normbvn[5,1] <- mean(bvn_2_1$model[,5]);normbvn[5,2] <- sd(bvn_2_1$model[,5])
normbvn[6,1] <- mean(bvn_4_1$model[,5]);normbvn[6,2] <- sd(bvn_4_1$model[,5])
normbvn[7,1] <- mean(bvn_2_1$model[,6]);normbvn[7,2] <- sd(bvn_2_1$model[,6])
normbvn[8,1] <- mean(bvn_4_1$model[,6]);normbvn[8,2] <- sd(bvn_4_1$model[,6])
normbvn[9,1] <- mean(bvn_2_1$model[,10]);normbvn[9,2] <- sd(bvn_2_1$model[,10])
normbvn[10,1] <- mean(bvn_4_1$model[,10]);normbvn[10,2] <- sd(bvn_4_1$model[,10])
normbvn[11,1] <- mean(bvn_2_1$model[,11]);normbvn[11,2] <- sd(bvn_2_1$model[,11])
normbvn[12,1] <- mean(bvn_4_1$model[,11]);normbvn[12,2] <- sd(bvn_4_1$model[,11])
normbvn[13,1] <- mean(bvn_2_1$model[,12]);normbvn[13,2] <- sd(bvn_2_1$model[,12])
normbvn[14,1] <- mean(bvn_4_1$model[,12]);normbvn[14,2] <- sd(bvn_4_1$model[,12])
normbvn[15,1] <- mean(bvn_2_1$model[,13]);normbvn[15,2] <- sd(bvn_2_1$model[,13]) 
normbvn[16,1] <- mean(bvn_4_1$model[,13]);normbvn[16,2] <- sd(bvn_4_1$model[,13])
normbvn[17,1] <- mean(bvn_2_1$model[,14]);normbvn[17,2] <- sd(bvn_2_1$model[,14])
normbvn[18,1] <- mean(bvn_4_1$model[,14]);normbvn[18,2] <- sd(bvn_4_1$model[,14])
normbvn[19,1] <- mean(bvn_2_1$model[,15]);normbvn[19,2] <- sd(bvn_2_1$model[,15])
normbvn[20,1] <- mean(bvn_4_1$model[,15]);normbvn[20,2] <- sd(bvn_4_1$model[,15])


#normdpmm[1,4] <- 0
normdpmm[1,1] <- mean(dpmm_4_1$model[,1]);normdpmm[1,2] <- sd(dpmm_4_1$model[,1])
#normdpmm[3,4] <- 0
normdpmm[2,1] <- mean(dpmm_4_1$model[,2]);normdpmm[2,2] <- sd(dpmm_4_1$model[,2])
#normdpmm[5,4] <- 0
normdpmm[3,1] <- mean(dpmm_4_1$model[,3]);normdpmm[3,2] <- sd(dpmm_4_1$model[,3])
#normdpmm[7,4] <- 0
normdpmm[4,1] <- mean(dpmm_4_1$model[,4]);normdpmm[4,2] <- sd(dpmm_4_1$model[,4])
#normdpmm[9,4] <- 0
normdpmm[5,1] <- mean(dpmm_4_1$model[,5]);normdpmm[5,2] <- sd(dpmm_4_1$model[,5])
#normdpmm[11,4] <- 0
normdpmm[6,1] <- mean(dpmm_4_1$model[,6]);normdpmm[6,2] <- sd(dpmm_4_1$model[,6])
#normdpmm[13,4] <- 0
normdpmm[7,1] <- mean(dpmm_4_1$model[,7]);normdpmm[7,2] <- sd(dpmm_4_1$model[,7])
#normdpmm[15,4] <- 0
normdpmm[8,1] <- mean(dpmm_4_1$model[,8]);normdpmm[8,2] <- sd(dpmm_4_1$model[,8])
#normdpmm[17,4] <- 0
normdpmm[9,1] <- mean(dpmm_4_1$model[,9]);normdpmm[9,2] <- sd(dpmm_4_1$model[,9])
#normdpmm[19,4] <- 0
normdpmm[10,1] <- mean(dpmm_4_1$model[,10]);normdpmm[10,2] <- sd(dpmm_4_1$model[,10])

#truth
norm0[,4] <- true20
norm1[,4] <- true2
normbvn[,4] <- true2
normdpmm[,4] <- true2d

#bias
norm0[,3] <- norm0[,1]-true20
norm1[,3] <- norm1[,1]-true2
normbvn[,3] <- normbvn[,1]-true2
normdpmm[,3] <- normdpmm[,1]-true2d
#--------------------------------------------------------

#-----------------------------------------------
#summary table for parameters skew errors

skew0 <- matrix(0,nrow=16,ncol=4)
skew1 <- matrix(0,nrow=20,ncol=4)
skewbvn <- matrix(0,nrow=20,ncol=4)
skewdpmm <- matrix(0,nrow=10,ncol=4)


skew0[1,1] <- mean(simple_2_2$model0[,11]);skew0[1,2] <- sd(simple_2_2$model0[,11])
skew0[2,1] <- mean(simple_4_2$model0[,11]);skew0[2,2] <- sd(simple_4_2$model0[,11])
skew0[3,1] <- mean(simple_2_2$model0[,12]);skew0[3,2] <- sd(simple_2_2$model0[,12])
skew0[4,1] <- mean(simple_4_2$model0[,12]);skew0[4,2] <- sd(simple_4_2$model0[,12])
skew0[9,1] <- mean(simple_2_2$model0[,5]);skew0[9,2] <- sd(simple_2_2$model0[,5])
skew0[10,1] <- mean(simple_4_2$model0[,5]);skew0[10,2] <- sd(simple_4_2$model0[,5])
skew0[7,1] <- mean(simple_2_2$model0[,6]);skew0[7,2] <- sd(simple_2_2$model0[,6])
skew0[8,1] <- mean(simple_4_2$model0[,6]);skew0[8,2] <- sd(simple_4_2$model0[,6])
skew0[5,1] <- mean(simple_2_2$model0[,7]);skew0[5,2] <- sd(simple_2_2$model0[,7])
skew0[6,1] <- mean(simple_4_2$model0[,7]);skew0[6,2] <- sd(simple_4_2$model0[,7])
skew0[15,1] <- mean(simple_2_2$model0[,8]);skew0[15,2] <- sd(simple_2_2$model0[,8])  
skew0[16,1] <- mean(simple_4_2$model0[,8]);skew0[16,2] <- sd(simple_4_2$model0[,8])
skew0[13,1] <- mean(simple_2_2$model0[,9]);skew0[13,2] <- sd(simple_2_2$model0[,9])
skew0[14,1] <- mean(simple_4_2$model0[,9]);skew0[14,2] <- sd(simple_4_2$model0[,9])
skew0[11,1] <- mean(simple_2_2$model0[,10]);skew0[11,2] <- sd(simple_2_2$model0[,10])
skew0[12,1] <- mean(simple_4_2$model0[,10]);skew0[12,2] <- sd(simple_4_2$model0[,10])


skew1[1,1] <- mean(simple_2_2$model1[,7]);skew1[1,2] <- sd(simple_2_2$model1[,7])
skew1[2,1] <- mean(simple_4_2$model1[,7]);skew1[2,2] <- sd(simple_4_2$model1[,7])
skew1[3,1] <- mean(simple_2_2$model1[,8]);skew1[3,2] <- sd(simple_2_2$model1[,8])
skew1[4,1] <- mean(simple_4_2$model1[,8]);skew1[4,2] <- sd(simple_4_2$model1[,8])
skew1[5,1] <- mean(simple_2_2$model1[,9]);skew1[5,2] <- sd(simple_2_2$model1[,9])
skew1[6,1] <- mean(simple_4_2$model1[,9]);skew1[6,2] <- sd(simple_4_2$model1[,9])
skew1[7,1] <- mean(simple_2_2$model1[,10]);skew1[7,2] <- sd(simple_2_2$model1[,10])
skew1[8,1] <- mean(simple_4_2$model1[,10]);skew1[8,2] <- sd(simple_4_2$model1[,10])
skew1[9,1] <- mean(simple_2_2$model1[,14]);skew1[9,2] <- sd(simple_2_2$model1[,14])
skew1[10,1] <- mean(simple_4_2$model1[,14]);skew1[10,2] <- sd(simple_4_2$model1[,14])
skew1[11,1] <- mean(simple_2_2$model1[,15]);skew1[11,2] <- sd(simple_2_2$model1[,15])
skew1[12,1] <- mean(simple_4_2$model1[,15]);skew1[12,2] <- sd(simple_4_2$model1[,15])
skew1[13,1] <- mean(simple_2_2$model1[,16]);skew1[13,2] <- sd(simple_2_2$model1[,16])
skew1[14,1] <- mean(simple_4_2$model1[,16]);skew1[14,2] <- sd(simple_4_2$model1[,16])
skew1[15,1] <- mean(simple_2_2$model1[,17]);skew1[15,2] <- sd(simple_2_2$model1[,17])  
skew1[16,1] <- mean(simple_4_2$model1[,17]);skew1[16,2] <- sd(simple_4_2$model1[,17])
skew1[17,1] <- mean(simple_2_2$model1[,18]);skew1[17,2] <- sd(simple_2_2$model1[,18])
skew1[18,1] <- mean(simple_4_2$model1[,18]);skew1[18,2] <- sd(simple_4_2$model1[,18])
skew1[19,1] <- mean(simple_2_2$model1[,19]);skew1[19,2] <- sd(simple_2_2$model1[,19])
skew1[20,1] <- mean(simple_4_2$model1[,19]);skew1[20,2] <- sd(simple_4_2$model1[,19])

skewbvn[1,1] <- mean(bvn_2_2$model[,3]);skewbvn[1,2] <- sd(bvn_2_2$model[,3])
skewbvn[2,1] <- mean(bvn_4_2$model[,3]);skewbvn[2,2] <- sd(bvn_4_2$model[,3])
skewbvn[3,1] <- mean(bvn_2_2$model[,4]);skewbvn[3,2] <- sd(bvn_2_2$model[,4])
skewbvn[4,1] <- mean(bvn_4_2$model[,4]);skewbvn[4,2] <- sd(bvn_4_2$model[,4])
skewbvn[5,1] <- mean(bvn_2_2$model[,5]);skewbvn[5,2] <- sd(bvn_2_2$model[,5])
skewbvn[6,1] <- mean(bvn_4_2$model[,5]);skewbvn[6,2] <- sd(bvn_4_2$model[,5])
skewbvn[7,1] <- mean(bvn_2_2$model[,6]);skewbvn[7,2] <- sd(bvn_2_2$model[,6])
skewbvn[8,1] <- mean(bvn_4_2$model[,6]);skewbvn[8,2] <- sd(bvn_4_2$model[,6])
skewbvn[9,1] <- mean(bvn_2_2$model[,10]);skewbvn[9,2] <- sd(bvn_2_2$model[,10])
skewbvn[10,1] <- mean(bvn_4_2$model[,10]);skewbvn[10,2] <- sd(bvn_4_2$model[,10])
skewbvn[11,1] <- mean(bvn_2_2$model[,11]);skewbvn[11,2] <- sd(bvn_2_2$model[,11])
skewbvn[12,1] <- mean(bvn_4_2$model[,11]);skewbvn[12,2] <- sd(bvn_4_2$model[,11])
skewbvn[13,1] <- mean(bvn_2_2$model[,12]);skewbvn[13,2] <- sd(bvn_2_2$model[,12])
skewbvn[14,1] <- mean(bvn_4_2$model[,12]);skewbvn[14,2] <- sd(bvn_4_2$model[,12])
skewbvn[15,1] <- mean(bvn_2_2$model[,13]);skewbvn[15,2] <- sd(bvn_2_2$model[,13]) 
skewbvn[16,1] <- mean(bvn_4_2$model[,13]);skewbvn[16,2] <- sd(bvn_4_2$model[,13])
skewbvn[17,1] <- mean(bvn_2_2$model[,14]);skewbvn[17,2] <- sd(bvn_2_2$model[,14])
skewbvn[18,1] <- mean(bvn_4_2$model[,14]);skewbvn[18,2] <- sd(bvn_4_2$model[,14])
skewbvn[19,1] <- mean(bvn_2_2$model[,15]);skewbvn[19,2] <- sd(bvn_2_2$model[,15])
skewbvn[20,1] <- mean(bvn_4_2$model[,15]);skewbvn[20,2] <- sd(bvn_4_2$model[,15])


#skewdpmm[1,4] <- 0
skewdpmm[1,1] <- mean(dpmm_4_2$model[,1]);skewdpmm[1,2] <- sd(dpmm_4_1$model[,1])
#skewdpmm[3,4] <- 0
skewdpmm[2,1] <- mean(dpmm_4_2$model[,2]);skewdpmm[2,2] <- sd(dpmm_4_2$model[,2])
#skewdpmm[5,4] <- 0
skewdpmm[3,1] <- mean(dpmm_4_2$model[,3]);skewdpmm[3,2] <- sd(dpmm_4_2$model[,3])
#skewdpmm[7,4] <- 0
skewdpmm[4,1] <- mean(dpmm_4_2$model[,4]);skewdpmm[4,2] <- sd(dpmm_4_2$model[,4])
#skewdpmm[9,4] <- 0
skewdpmm[5,1] <- mean(dpmm_4_2$model[,5]);skewdpmm[5,2] <- sd(dpmm_4_2$model[,5])
#skewdpmm[11,4] <- 0
skewdpmm[6,1] <- mean(dpmm_4_2$model[,6]);skewdpmm[6,2] <- sd(dpmm_4_2$model[,6])
#skewdpmm[13,4] <- 0
skewdpmm[7,1] <- mean(dpmm_4_2$model[,7]);skewdpmm[7,2] <- sd(dpmm_4_2$model[,7])
#skewdpmm[15,4] <- 0
skewdpmm[8,1] <- mean(dpmm_4_2$model[,8]);skewdpmm[8,2] <- sd(dpmm_4_2$model[,8])
#skewdpmm[17,4] <- 0
skewdpmm[9,1] <- mean(dpmm_4_2$model[,9]);skewdpmm[9,2] <- sd(dpmm_4_2$model[,9])
#skewdpmm[19,4] <- 0
skewdpmm[10,1] <- mean(dpmm_4_2$model[,10]);skewdpmm[10,2] <- sd(dpmm_4_2$model[,10])

#truth
skew0[,4] <- true20
skew1[,4] <- true2
skewbvn[,4] <- true2
skewdpmm[,4] <- true2d

#bias
skew0[,3] <- skew0[,1]-true20
skew1[,3] <- skew1[,1]-true2
skewbvn[,3] <- skewbvn[,1]-true2
skewdpmm[,3] <- skewdpmm[,1]-true2d


#-------------------------------------------


#----------------------------------------#
#same but for bimodal data

bimodal0 <- matrix(0,nrow=16,ncol=4)
bimodal1 <- matrix(0,nrow=20,ncol=4)
bimodalbvn <- matrix(0,nrow=20,ncol=4)
bimodaldpmm <- matrix(0,nrow=10,ncol=4)


bimodal0[1,1] <- mean(simple_2_3$model0[,11]);bimodal0[1,2] <- sd(simple_2_3$model0[,11])
bimodal0[2,1] <- mean(simple_4_3$model0[,11]);bimodal0[2,2] <- sd(simple_4_3$model0[,11])
bimodal0[3,1] <- mean(simple_2_3$model0[,12]);bimodal0[3,2] <- sd(simple_2_3$model0[,12])
bimodal0[4,1] <- mean(simple_4_3$model0[,12]);bimodal0[4,2] <- sd(simple_4_3$model0[,12])
bimodal0[9,1] <- mean(simple_2_3$model0[,5]);bimodal0[9,2] <- sd(simple_2_3$model0[,5])
bimodal0[10,1] <- mean(simple_4_3$model0[,5]);bimodal0[10,2] <- sd(simple_4_3$model0[,5])
bimodal0[7,1] <- mean(simple_2_3$model0[,6]);bimodal0[7,2] <- sd(simple_2_3$model0[,6])
bimodal0[8,1] <- mean(simple_4_3$model0[,6]);bimodal0[8,2] <- sd(simple_4_3$model0[,6])
bimodal0[5,1] <- mean(simple_2_3$model0[,7]);bimodal0[5,2] <- sd(simple_2_3$model0[,7])
bimodal0[6,1] <- mean(simple_4_3$model0[,7]);bimodal0[6,2] <- sd(simple_4_3$model0[,7])
bimodal0[15,1] <- mean(simple_2_3$model0[,8]);bimodal0[15,2] <- sd(simple_2_3$model0[,8])  
bimodal0[16,1] <- mean(simple_4_3$model0[,8]);bimodal0[16,2] <- sd(simple_4_3$model0[,8])
bimodal0[13,1] <- mean(simple_2_3$model0[,9]);bimodal0[13,2] <- sd(simple_2_3$model0[,9])
bimodal0[14,1] <- mean(simple_4_3$model0[,9]);bimodal0[14,2] <- sd(simple_4_3$model0[,9])
bimodal0[11,1] <- mean(simple_2_3$model0[,10]);bimodal0[11,2] <- sd(simple_2_3$model0[,10])
bimodal0[12,1] <- mean(simple_4_3$model0[,10]);bimodal0[12,2] <- sd(simple_4_3$model0[,10])


bimodal1[1,1] <- mean(simple_2_3$model1[,7]);bimodal1[1,2] <- sd(simple_2_3$model1[,7])
bimodal1[2,1] <- mean(simple_4_3$model1[,7]);bimodal1[2,2] <- sd(simple_4_3$model1[,7])
bimodal1[3,1] <- mean(simple_2_3$model1[,8]);bimodal1[3,2] <- sd(simple_2_3$model1[,8])
bimodal1[4,1] <- mean(simple_4_3$model1[,8]);bimodal1[4,2] <- sd(simple_4_3$model1[,8])
bimodal1[5,1] <- mean(simple_2_3$model1[,9]);bimodal1[5,2] <- sd(simple_2_3$model1[,9])
bimodal1[6,1] <- mean(simple_4_3$model1[,9]);bimodal1[6,2] <- sd(simple_4_3$model1[,9])
bimodal1[7,1] <- mean(simple_2_3$model1[,10]);bimodal1[7,2] <- sd(simple_2_3$model1[,10])
bimodal1[8,1] <- mean(simple_4_3$model1[,10]);bimodal1[8,2] <- sd(simple_4_3$model1[,10])
bimodal1[9,1] <- mean(simple_2_3$model1[,14]);bimodal1[9,2] <- sd(simple_2_3$model1[,14])
bimodal1[10,1] <- mean(simple_4_3$model1[,14]);bimodal1[10,2] <- sd(simple_4_3$model1[,14])
bimodal1[11,1] <- mean(simple_2_3$model1[,15]);bimodal1[11,2] <- sd(simple_2_3$model1[,15])
bimodal1[12,1] <- mean(simple_4_3$model1[,15]);bimodal1[12,2] <- sd(simple_4_3$model1[,15])
bimodal1[13,1] <- mean(simple_2_3$model1[,16]);bimodal1[13,2] <- sd(simple_2_3$model1[,16])
bimodal1[14,1] <- mean(simple_4_3$model1[,16]);bimodal1[14,2] <- sd(simple_4_3$model1[,16])
bimodal1[15,1] <- mean(simple_2_3$model1[,17]);bimodal1[15,2] <- sd(simple_2_3$model1[,17])  
bimodal1[16,1] <- mean(simple_4_3$model1[,17]);bimodal1[16,2] <- sd(simple_4_3$model1[,17])
bimodal1[17,1] <- mean(simple_2_3$model1[,18]);bimodal1[17,2] <- sd(simple_2_3$model1[,18])
bimodal1[18,1] <- mean(simple_4_3$model1[,18]);bimodal1[18,2] <- sd(simple_4_3$model1[,18])
bimodal1[19,1] <- mean(simple_2_3$model1[,19]);bimodal1[19,2] <- sd(simple_2_3$model1[,19])
bimodal1[20,1] <- mean(simple_4_3$model1[,19]);bimodal1[20,2] <- sd(simple_4_3$model1[,19])

bimodalbvn[1,1] <- mean(bvn_2_3$model[,3]);bimodalbvn[1,2] <- sd(bvn_2_3$model[,3])
bimodalbvn[2,1] <- mean(bvn_4_3$model[,3]);bimodalbvn[2,2] <- sd(bvn_4_3$model[,3])
bimodalbvn[3,1] <- mean(bvn_2_3$model[,4]);bimodalbvn[3,2] <- sd(bvn_2_3$model[,4])
bimodalbvn[4,1] <- mean(bvn_4_3$model[,4]);bimodalbvn[4,2] <- sd(bvn_4_3$model[,4])
bimodalbvn[5,1] <- mean(bvn_2_3$model[,5]);bimodalbvn[5,2] <- sd(bvn_2_3$model[,5])
bimodalbvn[6,1] <- mean(bvn_4_3$model[,5]);bimodalbvn[6,2] <- sd(bvn_4_3$model[,5])
bimodalbvn[7,1] <- mean(bvn_2_3$model[,6]);bimodalbvn[7,2] <- sd(bvn_2_3$model[,6])
bimodalbvn[8,1] <- mean(bvn_4_3$model[,6]);bimodalbvn[8,2] <- sd(bvn_4_3$model[,6])
bimodalbvn[9,1] <- mean(bvn_2_3$model[,10]);bimodalbvn[9,2] <- sd(bvn_2_3$model[,10])
bimodalbvn[10,1] <- mean(bvn_4_3$model[,10]);bimodalbvn[10,2] <- sd(bvn_4_3$model[,10])
bimodalbvn[11,1] <- mean(bvn_2_3$model[,11]);bimodalbvn[11,2] <- sd(bvn_2_3$model[,11])
bimodalbvn[12,1] <- mean(bvn_4_3$model[,11]);bimodalbvn[12,2] <- sd(bvn_4_3$model[,11])
bimodalbvn[13,1] <- mean(bvn_2_3$model[,12]);bimodalbvn[13,2] <- sd(bvn_2_3$model[,12])
bimodalbvn[14,1] <- mean(bvn_4_3$model[,12]);bimodalbvn[14,2] <- sd(bvn_4_3$model[,12])
bimodalbvn[15,1] <- mean(bvn_2_3$model[,13]);bimodalbvn[15,2] <- sd(bvn_2_3$model[,13]) 
bimodalbvn[16,1] <- mean(bvn_4_3$model[,13]);bimodalbvn[16,2] <- sd(bvn_4_3$model[,13])
bimodalbvn[17,1] <- mean(bvn_2_3$model[,14]);bimodalbvn[17,2] <- sd(bvn_2_3$model[,14])
bimodalbvn[18,1] <- mean(bvn_4_3$model[,14]);bimodalbvn[18,2] <- sd(bvn_4_3$model[,14])
bimodalbvn[19,1] <- mean(bvn_2_3$model[,15]);bimodalbvn[19,2] <- sd(bvn_2_3$model[,15])
bimodalbvn[20,1] <- mean(bvn_4_3$model[,15]);bimodalbvn[20,2] <- sd(bvn_4_3$model[,15])


#bimodaldpmm[1,4] <- 0
bimodaldpmm[1,1] <- mean(dpmm_4_3$model[,1]);bimodaldpmm[1,2] <- sd(dpmm_4_3$model[,1])
#bimodaldpmm[3,4] <- 0
bimodaldpmm[2,1] <- mean(dpmm_4_3$model[,2]);bimodaldpmm[2,2] <- sd(dpmm_4_3$model[,2])
#bimodaldpmm[5,4] <- 0
bimodaldpmm[3,1] <- mean(dpmm_4_3$model[,3]);bimodaldpmm[3,2] <- sd(dpmm_4_3$model[,3])
#bimodaldpmm[7,4] <- 0
bimodaldpmm[4,1] <- mean(dpmm_4_3$model[,4]);bimodaldpmm[4,2] <- sd(dpmm_4_3$model[,4])
#bimodaldpmm[9,4] <- 0
bimodaldpmm[5,1] <- mean(dpmm_4_3$model[,5]);bimodaldpmm[5,2] <- sd(dpmm_4_3$model[,5])
#bimodaldpmm[11,4] <- 0
bimodaldpmm[6,1] <- mean(dpmm_4_3$model[,6]);bimodaldpmm[6,2] <- sd(dpmm_4_3$model[,6])
#bimodaldpmm[13,4] <- 0
bimodaldpmm[7,1] <- mean(dpmm_4_3$model[,7]);bimodaldpmm[7,2] <- sd(dpmm_4_3$model[,7])
#bimodaldpmm[15,4] <- 0
bimodaldpmm[8,1] <- mean(dpmm_4_3$model[,8]);bimodaldpmm[8,2] <- sd(dpmm_4_3$model[,8])
#bimodaldpmm[17,4] <- 0
bimodaldpmm[9,1] <- mean(dpmm_4_3$model[,9]);bimodaldpmm[9,2] <- sd(dpmm_4_3$model[,9])
#bimodaldpmm[19,4] <- 0
bimodaldpmm[10,1] <- mean(dpmm_4_3$model[,10]);bimodaldpmm[10,2] <- sd(dpmm_4_3$model[,10])

#truth
bimodal0[,4] <- true20
bimodal1[,4] <- true2
bimodalbvn[,4] <- true2
bimodaldpmm[,4] <- true2d

#bias
bimodal0[,3] <- bimodal0[,1]-true20
bimodal1[,3] <- bimodal1[,1]-true2
bimodalbvn[,3] <- bimodalbvn[,1]-true2
bimodaldpmm[,3] <- bimodaldpmm[,1]-true2d

#-------------------------------------

modelnames <- c("", "Replicates", "Naive","Linear","BVN","DPMM","Truth")
tablenames <- c("", "Replicates","Mean Est","Std Err","Bias","Truth")
errors <- c("Normal", "", "Skewed", "", "Bimodal", "")
parameters <- c("$\\sigma_{yee}$","","$\\sigma_{yes}$","","$\\sigma_{wee}$","","$\\sigma_{wes}$","",
                "$\\gamma_{1,ee}$","", "$\\gamma_{2,ee}$","","$\\gamma_{3,ee}$","","$\\gamma_{1,es}$","",
                "$\\gamma_{2,es}$","","$\\gamma_{3,es}$","")
parameters0 <- c("$\\sigma_{yee}$","","$\\sigma_{yes}$","",
                "$\\gamma_{1,ee}$","", "$\\gamma_{2,ee}$","","$\\gamma_{3,ee}$","","$\\gamma_{1,es}$","",
                "$\\gamma_{2,es}$","","$\\gamma_{3,es}$","")
parametersd <- c("$\\sigma_{yee}$","$\\sigma_{yes}$","$\\sigma_{wee}$","$\\sigma_{wes}$",
                "$\\gamma_{1,ee}$", "$\\gamma_{2,ee}$","$\\gamma_{3,ee}$","$\\gamma_{1,es}$",
                "$\\gamma_{2,es}$","$\\gamma_{3,es}$")

c1_6 <- rep(c(2,4),3)
c1_20 <- rep(c(2,4),10)
c1_10 <- rep(c(4),10)
c1_16 <- rep(c(2,4),8)


pmateedf <- data.frame(cbind(errors,c1_6,round(pmatee,digits=2)))
pmatesdf <- data.frame(cbind(errors,c1_6,round(pmates,digits=2)))

names(pmateedf) <- modelnames[1:6]
names(pmatesdf) <- modelnames[1:6]

norm0df <- data.frame(t(cbind(c1_16,round(norm0,digits=2))))
norm1df <- data.frame(t(cbind(c1_20,round(norm1,digits=2))))
normbvndf <- data.frame(t(cbind(c1_20,round(normbvn,digits=2))))
normdpmmdf <- data.frame(t(cbind(c1_10,round(normdpmm,digits=2))))

skew0df <- data.frame(t(cbind(c1_16,round(skew0,digits=2))))
skew1df <- data.frame(t(cbind(c1_20,round(skew1,digits=2))))
skewbvndf <- data.frame(t(cbind(c1_20,round(skewbvn,digits=2))))
skewdpmmdf <- data.frame(t(cbind(c1_10,round(skewdpmm,digits=2))))

bimodal0df <- data.frame(t(cbind(c1_16,round(bimodal0,digits=2))))
bimodal1df <- data.frame(t(cbind(c1_20,round(bimodal1,digits=2))))
bimodalbvndf <- data.frame(t(cbind(c1_20,round(bimodalbvn,digits=2))))
bimodaldpmmdf <- data.frame(t(cbind(c1_10,round(bimodaldpmm,digits=2))))



names(norm0df) <- parameters0 
names(norm1df) <- parameters 
names(normbvndf) <- parameters 
names(normdpmmdf) <- parametersd 

names(skew0df) <- parameters0
names(skew1df) <- parameters 
names(skewbvndf) <- parameters 
names(skewdpmmdf) <- parametersd 

names(bimodal0df) <- parameters0 
names(bimodal1df) <- parameters 
names(bimodalbvndf) <- parameters 
names(bimodaldpmmdf) <- parametersd 

rownames(norm0df) <- tablenames[2:6] 
rownames(norm1df) <- tablenames[2:6]  
rownames(normbvndf) <- tablenames[2:6]  
rownames(normdpmmdf) <- tablenames[2:6]  

rownames(skew0df) <- tablenames[2:6] 
rownames(skew1df) <- tablenames[2:6]  
rownames(skewbvndf) <- tablenames[2:6]  
rownames(skewdpmmdf) <- tablenames[2:6]  

rownames(bimodal0df) <- tablenames[2:6]  
rownames(bimodal1df) <- tablenames[2:6]  
rownames(bimodalbvndf) <- tablenames[2:6]  
rownames(bimodaldpmmdf) <- tablenames[2:6]  

mdat0 <- matrix(c(rep(0,17),rep(2,17*4)),ncol=17,byrow=T)
mdat <- matrix(c(rep(0,21),rep(2,21*4)),ncol=21,byrow=T)
mdatd <- matrix(c(rep(0,11),rep(2,11*4)),ncol=11,byrow=T)


print(xtable(pmateedf,caption="PMSE for EE Regression",label="pmseee",align="lll|llll"),sanitize.text.function = function(x){x},include.rownames=FALSE)
print(xtable(pmatesdf,caption="PMSE for $\\Delta$ES Regression",label="pmsees",align="lll|llll"),sanitize.text.function = function(x){x},include.rownames=FALSE)

print(xtable((norm0df),digits=mdat0,align="r|ll|ll|ll|ll|ll|ll|ll|ll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(norm0df)),floating=FALSE)
print(xtable((norm1df),digits=mdat,align="r|ll|ll|ll|ll|ll|ll|ll|ll|ll|ll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(norm1df)),floating=FALSE)
print(xtable((normbvndf),digits=mdat,align="r|ll|ll|ll|ll|ll|ll|ll|ll|ll|ll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(normbvndf)),floating=FALSE)
print(xtable((normdpmmdf),digits=mdatd,align="r|llllllllll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(normdpmmdf)),floating=FALSE)

print(xtable((skew0df),digits=mdat0,align="r|ll|ll|ll|ll|ll|ll|ll|ll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(skew0df)),floating=FALSE)
print(xtable((skew1df),digits=mdat,align="r|ll|ll|ll|ll|ll|ll|ll|ll|ll|ll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(skew1df)),floating=FALSE)
print(xtable((skewbvndf),digits=mdat,align="r|ll|ll|ll|ll|ll|ll|ll|ll|ll|ll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(skewbvndf)),floating=FALSE)
print(xtable((skewdpmmdf),digits=mdatd,align="r|llllllllll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(skewdpmmdf)),floating=FALSE)

print(xtable((bimodal0df),digits=mdat0,align="r|ll|ll|ll|ll|ll|ll|ll|ll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(bimodal0df)),floating=FALSE)
print(xtable((bimodal1df),digits=mdat,align="r|ll|ll|ll|ll|ll|ll|ll|ll|ll|ll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(bimodal1df)),floating=FALSE)
print(xtable((bimodalbvndf),digits=mdat,align="r|ll|ll|ll|ll|ll|ll|ll|ll|ll|ll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(bimodalbvndf)),floating=FALSE)
print(xtable((bimodaldpmmdf),digits=mdatd,align="r|llllllllll"),sanitize.text.function = function(x){x},include.rownames=TRUE,hline.after = c(-1,0,1,nrow(bimodaldpmmdf)),floating=FALSE)

#\subcaption*{sdjf}
#\scalebox{0.8}{
#
