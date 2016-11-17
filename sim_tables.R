library(xtable)
load("C:/Users/dcries/bvn_2_1.RData")
load("C:/Users/dcries/bvn_2_2.RData")
load("C:/Users/dcries/bvn_2_3.RData")
load("C:/Users/dcries/bvn_4_1.RData")
load("C:/Users/dcries/bvn_4_2.RData")
load("C:/Users/dcries/bvn_4_3.RData")

load("C:/Users/dcries/dpmm_4_1.RData")
load("C:/Users/dcries/dpmm_4_2.RData")
load("C:/Users/dcries/dpmm_4_3.RData")


out <- bvn_4_1
pmse <- out$pmse
model <- out$model[,-c(1,2,7,8,9)]
indcheck <- out$indcheck


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


pmatee <- matrix(0,nrow=4,ncol=6)
pmatee[1,1] <- mean(bvn_2_1$pmse[,1])
pmatee[1,2] <- mean(bvn_4_1$pmse[,1])
pmatee[1,3] <- mean(bvn_2_2$pmse[,1])
pmatee[1,4] <- mean(bvn_4_3$pmse[,1])
pmatee[1,5] <- mean(bvn_2_3$pmse[,1])
pmatee[1,6] <- mean(bvn_4_2$pmse[,1])
pmatee[2,1] <- mean(bvn_2_1$pmse[,1])
pmatee[2,2] <- mean(bvn_4_1$pmse[,1])
pmatee[2,3] <- mean(bvn_2_2$pmse[,1])
pmatee[2,4] <- mean(bvn_4_2$pmse[,1])
pmatee[2,5] <- mean(bvn_2_3$pmse[,1])
pmatee[2,6] <- mean(bvn_4_3$pmse[,1])
pmatee[3,1] <- mean(bvn_2_1$pmse[,1])
pmatee[3,2] <- mean(bvn_4_1$pmse[,1])
pmatee[3,3] <- mean(bvn_2_2$pmse[,1])
pmatee[3,4] <- mean(bvn_4_2$pmse[,1])
pmatee[3,5] <- mean(bvn_2_3$pmse[,1])
pmatee[3,6] <- mean(bvn_4_3$pmse[,1])
pmatee[4,1] <- mean(bvn_2_1$pmse[,1])
pmatee[4,2] <- mean(bvn_4_1$pmse[,1])
pmatee[4,3] <- mean(bvn_2_2$pmse[,1])
pmatee[4,4] <- mean(bvn_4_2$pmse[,1])
pmatee[4,5] <- mean(bvn_2_3$pmse[,1])
pmatee[4,6] <- mean(bvn_4_3$pmse[,1])

pmates <- matrix(0,nrow=4,ncol=6)
pmates[1,1] <- mean(bvn_2_1$pmse[,2])
pmates[1,2] <- mean(bvn_4_1$pmse[,2])
pmates[1,3] <- mean(bvn_2_2$pmse[,2])
pmates[1,4] <- mean(bvn_4_3$pmse[,2])
pmates[1,5] <- mean(bvn_2_3$pmse[,2])
pmates[1,6] <- mean(bvn_4_2$pmse[,2])
pmates[2,1] <- mean(bvn_2_1$pmse[,2])
pmates[2,2] <- mean(bvn_4_1$pmse[,2])
pmates[2,3] <- mean(bvn_2_2$pmse[,2])
pmates[2,4] <- mean(bvn_4_2$pmse[,2])
pmates[2,5] <- mean(bvn_2_3$pmse[,2])
pmates[2,6] <- mean(bvn_4_3$pmse[,2])
pmates[3,1] <- mean(bvn_2_1$pmse[,2])
pmates[3,2] <- mean(bvn_4_1$pmse[,2])
pmates[3,3] <- mean(bvn_2_2$pmse[,2])
pmates[3,4] <- mean(bvn_4_2$pmse[,2])
pmates[3,5] <- mean(bvn_2_3$pmse[,2])
pmates[3,6] <- mean(bvn_4_3$pmse[,2])
pmates[4,1] <- mean(bvn_2_1$pmse[,2])
pmates[4,2] <- mean(bvn_4_1$pmse[,2])
pmates[4,3] <- mean(bvn_2_2$pmse[,2])
pmates[4,4] <- mean(bvn_4_2$pmse[,2])
pmates[4,5] <- mean(bvn_2_3$pmse[,2])
pmates[4,6] <- mean(bvn_4_3$pmse[,2])
