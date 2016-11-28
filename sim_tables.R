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

true2 <- c(420,0,344,0,250,0,72.862,0,300,0,14,0,-7,0,-200,0,8,0,-5,0)
#-----------------------------------------------
#summary table for parameters normal errors

norm <- matrix(0,nrow=20,ncol=5)

norm[1,1] <- mean(simple_2_1$model0[,11])
norm[2,1] <- mean(simple_4_1$model0[,11])
norm[3,1] <- mean(simple_2_1$model0[,12])
norm[4,1] <- mean(simple_4_1$model0[,12])
norm[5,1] <- 0
norm[6,1] <- 0
norm[7,1] <- 0
norm[8,1] <- 0
norm[9,1] <- mean(simple_2_1$model0[,5])
norm[10,1] <- mean(simple_4_1$model0[,5])
norm[11,1] <- mean(simple_2_1$model0[,6])
norm[12,1] <- mean(simple_4_1$model0[,6])
norm[13,1] <- mean(simple_2_1$model0[,7])
norm[14,1] <- mean(simple_4_1$model0[,7])
norm[15,1] <- mean(simple_2_1$model0[,8]) 
norm[16,1] <- mean(simple_4_1$model0[,8])
norm[17,1] <- mean(simple_2_1$model0[,9])
norm[18,1] <- mean(simple_4_1$model0[,9])
norm[19,1] <- mean(simple_2_1$model0[,10])
norm[20,1] <- mean(simple_4_1$model0[,10])


norm[1,2] <- mean(simple_2_1$model1[,7])
norm[2,2] <- mean(simple_4_1$model1[,7])
norm[3,2] <- mean(simple_2_1$model1[,8])
norm[4,2] <- mean(simple_4_1$model1[,8])
norm[5,2] <- mean(simple_2_1$model1[,9])
norm[6,2] <- mean(simple_4_1$model1[,9])
norm[7,2] <- mean(simple_2_1$model1[,10])
norm[8,2] <- mean(simple_4_1$model1[,10])
norm[9,2] <- mean(simple_2_1$model1[,14])
norm[10,2] <- mean(simple_4_1$model1[,14])
norm[11,2] <- mean(simple_2_1$model1[,15])
norm[12,2] <- mean(simple_4_1$model1[,15])
norm[13,2] <- mean(simple_2_1$model1[,16])
norm[14,2] <- mean(simple_4_1$model1[,16])
norm[15,2] <- mean(simple_2_1$model1[,17]) 
norm[16,2] <- mean(simple_4_1$model1[,17])
norm[17,2] <- mean(simple_2_1$model1[,18])
norm[18,2] <- mean(simple_4_1$model1[,18])
norm[19,2] <- mean(simple_2_1$model1[,19])
norm[20,2] <- mean(simple_4_1$model1[,19])

norm[1,3] <- mean(bvn_2_1$model[,3])
norm[2,3] <- mean(bvn_4_1$model[,3])
norm[3,3] <- mean(bvn_2_1$model[,4])
norm[4,3] <- mean(bvn_4_1$model[,4])
norm[5,3] <- mean(bvn_2_1$model[,5])
norm[6,3] <- mean(bvn_4_1$model[,5])
norm[7,3] <- mean(bvn_2_1$model[,6])
norm[8,3] <- mean(bvn_4_1$model[,6])
norm[9,3] <- mean(bvn_2_1$model[,10])
norm[10,3] <- mean(bvn_4_1$model[,10])
norm[11,3] <- mean(bvn_2_1$model[,11])
norm[12,3] <- mean(bvn_4_1$model[,11])
norm[13,3] <- mean(bvn_2_1$model[,12])
norm[14,3] <- mean(bvn_4_1$model[,12])
norm[15,3] <- mean(bvn_2_1$model[,13]) 
norm[16,3] <- mean(bvn_4_1$model[,13])
norm[17,3] <- mean(bvn_2_1$model[,14])
norm[18,3] <- mean(bvn_4_1$model[,14])
norm[19,3] <- mean(bvn_2_1$model[,15])
norm[20,3] <- mean(bvn_4_1$model[,15])


norm[1,4] <- 0
norm[2,4] <- mean(dpmm_4_1$model[,1])
norm[3,4] <- 0
norm[4,4] <- mean(dpmm_4_1$model[,2])
norm[5,4] <- 0
norm[6,4] <- mean(dpmm_4_1$model[,3])
norm[7,4] <- 0
norm[8,4] <- mean(dpmm_4_1$model[,4])
norm[9,4] <- 0
norm[10,4] <- mean(dpmm_4_1$model[,5])
norm[11,4] <- 0
norm[12,4] <- mean(dpmm_4_1$model[,6])
norm[13,4] <- 0
norm[14,4] <- mean(dpmm_4_1$model[,7])
norm[15,4] <- 0
norm[16,4] <- mean(dpmm_4_1$model[,8])
norm[17,4] <- 0
norm[18,4] <- mean(dpmm_4_1$model[,9])
norm[19,4] <- 0
norm[20,4] <- mean(dpmm_4_1$model[,10])

#truth
norm[,5] <- true2

#-------------------------------------------

#----------------------------------------#
#same but for skewed data

skew <- matrix(0,nrow=20,ncol=5)

skew[1,1] <- mean(simple_2_2$model0[,11])
skew[2,1] <- mean(simple_4_2$model0[,11])
skew[3,1] <- mean(simple_2_2$model0[,12])
skew[4,1] <- mean(simple_4_2$model0[,12])
skew[5,1] <- 0
skew[6,1] <- 0
skew[7,1] <- 0
skew[8,1] <- 0
skew[9,1] <- mean(simple_2_2$model0[,5])
skew[10,1] <- mean(simple_4_2$model0[,5])
skew[11,1] <- mean(simple_2_2$model0[,6])
skew[12,1] <- mean(simple_4_2$model0[,6])
skew[13,1] <- mean(simple_2_2$model0[,7])
skew[14,1] <- mean(simple_4_2$model0[,7])
skew[15,1] <- mean(simple_2_2$model0[,8]) 
skew[16,1] <- mean(simple_4_2$model0[,8])
skew[17,1] <- mean(simple_2_2$model0[,9])
skew[18,1] <- mean(simple_4_2$model0[,9])
skew[19,1] <- mean(simple_2_2$model0[,10])
skew[20,1] <- mean(simple_4_2$model0[,10])


skew[1,2] <- mean(simple_2_2$model1[,7])
skew[2,2] <- mean(simple_4_2$model1[,7])
skew[3,2] <- mean(simple_2_2$model1[,8])
skew[4,2] <- mean(simple_4_2$model1[,8])
skew[5,2] <- mean(simple_2_2$model1[,9])
skew[6,2] <- mean(simple_4_2$model1[,9])
skew[7,2] <- mean(simple_2_2$model1[,10])
skew[8,2] <- mean(simple_4_2$model1[,10])
skew[9,2] <- mean(simple_2_2$model1[,14])
skew[10,2] <- mean(simple_4_2$model1[,14])
skew[11,2] <- mean(simple_2_2$model1[,15])
skew[12,2] <- mean(simple_4_2$model1[,15])
skew[13,2] <- mean(simple_2_2$model1[,16])
skew[14,2] <- mean(simple_4_2$model1[,16])
skew[15,2] <- mean(simple_2_2$model1[,17]) 
skew[16,2] <- mean(simple_4_2$model1[,17])
skew[17,2] <- mean(simple_2_2$model1[,18])
skew[18,2] <- mean(simple_4_2$model1[,18])
skew[19,2] <- mean(simple_2_2$model1[,19])
skew[20,2] <- mean(simple_4_2$model1[,19])

skew[1,3] <- mean(bvn_2_2$model[,3])
skew[2,3] <- mean(bvn_4_2$model[,3])
skew[3,3] <- mean(bvn_2_2$model[,4])
skew[4,3] <- mean(bvn_4_2$model[,4])
skew[5,3] <- mean(bvn_2_2$model[,5])
skew[6,3] <- mean(bvn_4_2$model[,5])
skew[7,3] <- mean(bvn_2_2$model[,6])
skew[8,3] <- mean(bvn_4_2$model[,6])
skew[9,3] <- mean(bvn_2_2$model[,10])
skew[10,3] <- mean(bvn_4_2$model[,10])
skew[11,3] <- mean(bvn_2_2$model[,11])
skew[12,3] <- mean(bvn_4_2$model[,11])
skew[13,3] <- mean(bvn_2_2$model[,12])
skew[14,3] <- mean(bvn_4_2$model[,12])
skew[15,3] <- mean(bvn_2_2$model[,13]) 
skew[16,3] <- mean(bvn_4_2$model[,13])
skew[17,3] <- mean(bvn_2_2$model[,14])
skew[18,3] <- mean(bvn_4_2$model[,14])
skew[19,3] <- mean(bvn_2_2$model[,15])
skew[20,3] <- mean(bvn_4_2$model[,15])


skew[1,4] <- 0
skew[2,4] <- mean(dpmm_4_2$model[,1])
skew[3,4] <- 0
skew[4,4] <- mean(dpmm_4_2$model[,2])
skew[5,4] <- 0
skew[6,4] <- mean(dpmm_4_2$model[,3])
skew[7,4] <- 0
skew[8,4] <- mean(dpmm_4_2$model[,4])
skew[9,4] <- 0
skew[10,4] <- mean(dpmm_4_2$model[,5])
skew[11,4] <- 0
skew[12,4] <- mean(dpmm_4_2$model[,6])
skew[13,4] <- 0
skew[14,4] <- mean(dpmm_4_2$model[,7])
skew[15,4] <- 0
skew[16,4] <- mean(dpmm_4_2$model[,8])
skew[17,4] <- 0
skew[18,4] <- mean(dpmm_4_2$model[,9])
skew[19,4] <- 0
skew[20,4] <- mean(dpmm_4_2$model[,10])

#truth
skew[,5] <- true2

#--------------------------------------------------------

#----------------------------------------#
#same but for bimodal data

bimodal <- matrix(0,nrow=20,ncol=5)

bimodal[1,1] <- mean(simple_2_3$model0[,11])
bimodal[2,1] <- mean(simple_4_3$model0[,11])
bimodal[3,1] <- mean(simple_2_3$model0[,12])
bimodal[4,1] <- mean(simple_4_3$model0[,12])
bimodal[5,1] <- 0
bimodal[6,1] <- 0
bimodal[7,1] <- 0
bimodal[8,1] <- 0
bimodal[9,1] <- mean(simple_2_3$model0[,5])
bimodal[10,1] <- mean(simple_4_3$model0[,5])
bimodal[11,1] <- mean(simple_2_3$model0[,6])
bimodal[12,1] <- mean(simple_4_3$model0[,6])
bimodal[13,1] <- mean(simple_2_3$model0[,7])
bimodal[14,1] <- mean(simple_4_3$model0[,7])
bimodal[15,1] <- mean(simple_2_3$model0[,8]) 
bimodal[16,1] <- mean(simple_4_3$model0[,8])
bimodal[17,1] <- mean(simple_2_3$model0[,9])
bimodal[18,1] <- mean(simple_4_3$model0[,9])
bimodal[19,1] <- mean(simple_2_3$model0[,10])
bimodal[20,1] <- mean(simple_4_3$model0[,10])


bimodal[1,2] <- mean(simple_2_3$model1[,7])
bimodal[2,2] <- mean(simple_4_3$model1[,7])
bimodal[3,2] <- mean(simple_2_3$model1[,8])
bimodal[4,2] <- mean(simple_4_3$model1[,8])
bimodal[5,2] <- mean(simple_2_3$model1[,9])
bimodal[6,2] <- mean(simple_4_3$model1[,9])
bimodal[7,2] <- mean(simple_2_3$model1[,10])
bimodal[8,2] <- mean(simple_4_3$model1[,10])
bimodal[9,2] <- mean(simple_2_3$model1[,14])
bimodal[10,2] <- mean(simple_4_3$model1[,14])
bimodal[11,2] <- mean(simple_2_3$model1[,15])
bimodal[12,2] <- mean(simple_4_3$model1[,15])
bimodal[13,2] <- mean(simple_2_3$model1[,16])
bimodal[14,2] <- mean(simple_4_3$model1[,16])
bimodal[15,2] <- mean(simple_2_3$model1[,17]) 
bimodal[16,2] <- mean(simple_4_3$model1[,17])
bimodal[17,2] <- mean(simple_2_3$model1[,18])
bimodal[18,2] <- mean(simple_4_3$model1[,18])
bimodal[19,2] <- mean(simple_2_3$model1[,19])
bimodal[20,2] <- mean(simple_4_3$model1[,19])

bimodal[1,3] <- mean(bvn_2_3$model[,3])
bimodal[2,3] <- mean(bvn_4_3$model[,3])
bimodal[3,3] <- mean(bvn_2_3$model[,4])
bimodal[4,3] <- mean(bvn_4_3$model[,4])
bimodal[5,3] <- mean(bvn_2_3$model[,5])
bimodal[6,3] <- mean(bvn_4_3$model[,5])
bimodal[7,3] <- mean(bvn_2_3$model[,6])
bimodal[8,3] <- mean(bvn_4_3$model[,6])
bimodal[9,3] <- mean(bvn_2_3$model[,10])
bimodal[10,3] <- mean(bvn_4_3$model[,10])
bimodal[11,3] <- mean(bvn_2_3$model[,11])
bimodal[12,3] <- mean(bvn_4_3$model[,11])
bimodal[13,3] <- mean(bvn_2_3$model[,12])
bimodal[14,3] <- mean(bvn_4_3$model[,12])
bimodal[15,3] <- mean(bvn_2_3$model[,13]) 
bimodal[16,3] <- mean(bvn_4_3$model[,13])
bimodal[17,3] <- mean(bvn_2_3$model[,14])
bimodal[18,3] <- mean(bvn_4_3$model[,14])
bimodal[19,3] <- mean(bvn_2_3$model[,15])
bimodal[20,3] <- mean(bvn_4_3$model[,15])


bimodal[1,4] <- 0
bimodal[2,4] <- mean(dpmm_4_3$model[,1])
bimodal[3,4] <- 0
bimodal[4,4] <- mean(dpmm_4_3$model[,2])
bimodal[5,4] <- 0
bimodal[6,4] <- mean(dpmm_4_3$model[,3])
bimodal[7,4] <- 0
bimodal[8,4] <- mean(dpmm_4_3$model[,4])
bimodal[9,4] <- 0
bimodal[10,4] <- mean(dpmm_4_3$model[,5])
bimodal[11,4] <- 0
bimodal[12,4] <- mean(dpmm_4_3$model[,6])
bimodal[13,4] <- 0
bimodal[14,4] <- mean(dpmm_4_3$model[,7])
bimodal[15,4] <- 0
bimodal[16,4] <- mean(dpmm_4_3$model[,8])
bimodal[17,4] <- 0
bimodal[18,4] <- mean(dpmm_4_3$model[,9])
bimodal[19,4] <- 0
bimodal[20,4] <- mean(dpmm_4_3$model[,10])

#truth
bimodal[,5] <- true2

#-------------------------------------

modelnames <- c("", "Replicates", "Naive","Linear","BVN","DPMM","Truth")
errors <- c("Normal", "", "Skewed", "", "Bimodal", "")
parameters <- c("$\\sigma_{yee}$","","$\\sigma_{yes}$","","$\\sigma_{wee}$","","$\\sigma_{wes}$","",
                "$\\gamma_{1,ee}$","", "$\\gamma_{2,ee}$","","$\\gamma_{3,ee}$","","$\\gamma_{1,es}$","",
                "$\\gamma_{2,es}$","","$\\gamma_{3,es}$","")
c1_6 <- rep(c(2,4),3)
c1_20 <- rep(c(2,4),10)

pmateedf <- data.frame(cbind(errors,c1_6,round(pmatee,digits=2)))
pmatesdf <- data.frame(cbind(errors,c1_6,round(pmates,digits=2)))

names(pmateedf) <- modelnames[1:6]
names(pmatesdf) <- modelnames[1:6]

normdf <- data.frame(cbind(parameters,c1_20,round(norm,digits=2)))
skewdf <- data.frame(cbind(parameters,c1_20,round(skew,digits=2)))
bimodaldf <- data.frame(cbind(parameters,c1_20,round(bimodal,digits=2)))

names(normdf) <- modelnames 
names(skewdf) <- modelnames 
names(bimodaldf) <- modelnames 

print(xtable(pmateedf,caption="dsfg",label="dfsf",align="lll|llll"),sanitize.text.function = function(x){x},include.rownames=FALSE)
print(xtable(pmatesdf,caption="dsfg",label="dfsf",align="lll|llll"),sanitize.text.function = function(x){x},include.rownames=FALSE)
print(xtable(normdf,caption="dsfg",label="dfsf",align="lll|llll"),sanitize.text.function = function(x){x},include.rownames=FALSE)
print(xtable(skewdf,caption="dsfg",label="dfsf",align="lll|llll"),sanitize.text.function = function(x){x},include.rownames=FALSE)
print(xtable(bimodaldf,caption="dsfg",label="dfsf",align="lll|llll"),sanitize.text.function = function(x){x},include.rownames=FALSE)

