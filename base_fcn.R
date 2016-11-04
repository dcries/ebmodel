#Rcpp::sourceCpp("\\\\my.files.iastate.edu\\Users\\dcries\\Desktop\\research\\rprograms\\dp_calc.cpp")

inverse <- function (f, lower = -100, upper = 100) {
  #function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
  function (y) optimize((function (x) abs(f(x) - y)), lower = lower, upper = upper,maximum=FALSE)$minimum
}

callibrate <- function(y,z,latentx,knotlocations,betacoef,gammacoef,nk,ng,ndraws,min=0,max=5500){
  #given a value of y, we want to calculate the corresponding value of x
  # y is value to predict with
  # z is vector of observed error free covariates corresponding to individual of interest
  # latent x is a matrix of latent x variables for all individuals, each col is individual and each 
  # row represents another draw
  # knotlocations is matrix of knot locations
  indicies <- sample(1:nrow(gammacoef),ndraws)
  prediction <- rep(0,ndraws)
  if(!is.matrix(betacoef)){
    betacoef <- matrix(0,nrow=nrow(gammacoef),ncol=1)
  }
  ct <- 1
  for(i in indicies){
    basis_matrix <- bs(latentx[i,],knots=knotlocations[i,1:nk[i]],intercept=TRUE)
    coefs <- gammacoef[i,1:ng[i]]
    my_inverse <- inverse(function(x) {return(coefs%*%(as.matrix(predict(basis_matrix,x))[1,]))},min,max)
    prediction[ct] <- my_inverse(y- t(z)%*%betacoef[i,]) 
    ct <- ct+1
  }
  return(prediction)
}

density_conv <- function(data,x,chains=3){
  #data: matrix output from JAGS
  # col represent parameters, each col has data from
  #all chains ordered, first m values in col are chain1, second m values in col
  #are chain2, etc
  #assuming mixture of normals was used, and only parameters are mu,pi,sigma
  #in that order
  
  #x: is the location to get density at
  #chains: is the number of chains
  h <- ncol(data)/3 #assuming only pi, mu, sigma are parameters
  #niter <- nrow(data)/chains
  
  mu <- data[,1:h]
  pi <- data[,((h+1):(2*h))]
  sigma <- data[,(2*h+1):(3*h)]
  m <- nrow(mu)/chains #number of mcmc iterations
  
  #temp <- 0
  temp <- rep(0,nrow(data))
  for(j in 1:h){
    temp <- temp + pi[,j]*dnorm(x,mu[,j],sigma[,j]) #[,j] are mx1 vectors
  }
  
  den <- matrix(temp,nrow=m,ncol=chains,byrow=FALSE)
  den <- as.data.frame(den)
  den$m <- 1:m
  for(k in 1:chains){
    names(den)[k] <- paste0("chain",k)
  }
  mden <- melt(den,id.vars="m")
  ggplot(data=mden) + geom_line(aes(x=m,y=value,colour=variable),alpha=0.7)
  
}

plot_dp <- function(real,post,evals){
  #obs is the vector of observed data
  #post is a matrix with mcmc output from JAGS, provide only the components
  #will be using, ie mu, pi, sigma
  #evals is a vector of length 3, with min and max values, and number
  # of terms to evaluate density
  real <- as.data.frame(real)
  h <- ncol(post)/3 #h is the number of components in the model, divide by 3 
  #since i'm assuming i only have mu, pi, and sigma
  eval <- seq(evals[1],evals[2],length.out=evals[3])
  
  mu <- post[,1:h]
  pi <- post[,((h+1):(2*h))]
  sigma <- sqrt(post[,(2*h+1):(3*h)])
  
  #   final <- matrix(rep(0,nrow(post)*length(eval)),ncol=nrow(post),nrow=length(eval))
  #   for(i in 1:nrow(post)){
  #     temp <- rep(0,length(eval))
  #     for(j in 1:h){
  #       temp <- temp + pi[i,j]*dnorm(eval,mu[i,j],sigma[i,j])
  #     }
  #     final[,i] <- temp
  #   }
  final <- dp_postdensity(mu,pi,sigma,eval)
  
  
  ci <- as.data.frame(t(apply(final,1,quantile,probs=c(0.025,0.5,0.975))))
  names(ci) <- c("low","med","up")
  ci$eval <- eval
  
  ggplot() + geom_histogram(data=real,aes(x=real,..density..)) + geom_line(data=ci,aes(y=med,x=eval),colour="red",size=1) +
    geom_ribbon(data=ci,aes(x = eval, ymin = low, ymax = up), fill = "blue", alpha = 0.4)
  
}

ppred_var <- function(post,truth){
  #function to check truth to posterior predictive distribution of 
  #variance term
  #assume density is normal
  #need ggplot
  #post is a vector of posterior draws
  #truth is true std dev
  n <- length(post)
  pred <- rnorm(n,0,post)
  loc <- seq(min(pred),max(pred),by=1)
  eval <- dnorm(loc,0,truth)
  dat1 <- data.frame(pred=pred)
  dat2 <- data.frame(loc=loc,eval=eval)
  ggplot() + geom_histogram(data=dat1,aes(x=pred,..density..)) + geom_line(data=dat2,aes(x=loc,y=eval),colour="red")
}

bivariate.normal <- function(x, mu, Sigma) {
  exp(-.5*t(x-mu)%*%solve(Sigma)%*%(x-mu))/sqrt(2*pi*det(Sigma))
}


basis <- function(x,knots){
  #computes basis function X matrix for cubic polynomial spline
  #x is vector of x values sorted least to greatest, and knots is knot locations
  #x <- sort(x)
  n <- length(x)
  K <- length(knots)
  B <- matrix(0,n,4+length(knots)) #cubic with intercept
  B[,1] <- rep(1,n)
  B[,2] <- x
  B[,3] <- x^2
  B[,4] <- x^3
  
  for(i in 1:K){
    zeros <- sum(x <= knots[i])
    #B[(1:zeros),(i+4)] <- rep(0,zeros)
    B[((zeros+1):n),(i+4)] <- (x[x>knots[i]]-knots[i])^3
  }
  return(B)
}


B.basis <- function(x,knots)
{
  delta <- knots[2]-knots[1]
  n <- length(x)
  K <- length(knots)
  B <- matrix(0,n,K+1)
  for (jj in 1:(K-1))
  {
    act.inds <- (1:n)[(x>=knots[jj])&(x<=knots[jj+1])]
    act.x <- x[act.inds]
    resc.x <- (act.x-knots[jj])/(knots[jj+1]-knots[jj])
    
    B[act.inds,jj] <- (1/2)*(1-resc.x)^2
    B[act.inds,jj+1] <- -(resc.x^2)+resc.x+1/2
    B[act.inds,jj+2] <- (resc.x^2)/2
  }
  return(B)
}

#eecurve <- function(x,L=4000,x0=2200,k=0.002){
eecurve <- function(x,L=3800,x0=2000,k=0.0019){
  
  #L is maximum
  #x0 is midpoint 
  return(2*x-(L/(1+exp(-k*(x-x0)))))
}

#escurve <- function(x,L=1000,x0=0,k=0.04){ #1000, 0.04
escurve <- function(x,L=600,x0=0,k=0.04){ #1000, 0.04
  
  #return(x+400*cos(0.003*x)) # 500, 0.03
  return(L/(1+exp(-k*(x-x0)))-L/2 + x)
}

ppred8 <- function(data,true){
  #plots 10 histograms of draws from post pred dist for 10 different 
  #interesting statistics, and compares to the truth (true vector of 10)
  p1 <- qplot(x=data[,1],geom="histogram") + geom_vline(xintercept=true[1],colour="red") + theme_bw()
  p2 <- qplot(x=data[,2],geom="histogram") + geom_vline(xintercept=true[2],colour="red") + theme_bw()
  p3 <- qplot(x=data[,3],geom="histogram") + geom_vline(xintercept=true[3],colour="red") + theme_bw()
  p4 <- qplot(x=data[,4],geom="histogram") + geom_vline(xintercept=true[4],colour="red") + theme_bw()
  p5 <- qplot(x=data[,5],geom="histogram") + geom_vline(xintercept=true[5],colour="red") + theme_bw()
  p6 <- qplot(x=data[,6],geom="histogram") + geom_vline(xintercept=true[6],colour="red") + theme_bw()
  p7 <- qplot(x=data[,7],geom="histogram") + geom_vline(xintercept=true[7],colour="red") + theme_bw()
  p8 <- qplot(x=data[,8],geom="histogram") + geom_vline(xintercept=true[8],colour="red") + theme_bw()
  # p9 <- qplot(x=data[,9],geom="histogram") + geom_vline(xintercept=true[9],colour="red") + theme_bw()
  # p10 <- qplot(x=data[,10],geom="histogram") + geom_vline(xintercept=true[10],colour="red") + theme_bw()
  # p11 <- qplot(x=data[,11],geom="histogram") + geom_vline(xintercept=true[11],colour="red") + theme_bw()
  # p12 <- qplot(x=data[,12],geom="histogram") + geom_vline(xintercept=true[12],colour="red") + theme_bw()
  # p13 <- qplot(x=data[,13],geom="histogram") + geom_vline(xintercept=true[13],colour="red") + theme_bw()
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8)
  #do.call(grid.arrange,p)
}

ppred16 <- function(data,true,xlabels=as.character(1:16)){
  #plots 10 histograms of draws from post pred dist for 10 different 
  #interesting statistics, and compares to the truth (true vector of 10)
  p1 <- qplot(x=data[,1],geom="histogram") + geom_vline(xintercept=true[1],colour="red") + xlab(xlabels[1]) + theme_bw()
  p2 <- qplot(x=data[,2],geom="histogram") + geom_vline(xintercept=true[2],colour="red") + xlab(xlabels[2]) + theme_bw()
  p3 <- qplot(x=data[,3],geom="histogram") + geom_vline(xintercept=true[3],colour="red") + xlab(xlabels[3]) + theme_bw()
  p4 <- qplot(x=data[,4],geom="histogram") + geom_vline(xintercept=true[4],colour="red") + xlab(xlabels[4]) + theme_bw()
  p5 <- qplot(x=data[,5],geom="histogram") + geom_vline(xintercept=true[5],colour="red") + xlab(xlabels[5]) + theme_bw()
  p6 <- qplot(x=data[,6],geom="histogram") + geom_vline(xintercept=true[6],colour="red") + xlab(xlabels[6]) + theme_bw()
  p7 <- qplot(x=data[,7],geom="histogram") + geom_vline(xintercept=true[7],colour="red")  + xlab(xlabels[7]) + theme_bw()
  p8 <- qplot(x=data[,8],geom="histogram") + geom_vline(xintercept=true[8],colour="red") + xlab(xlabels[8]) + theme_bw()
  p9 <- qplot(x=data[,9],geom="histogram") + geom_vline(xintercept=true[9],colour="red")  + xlab(xlabels[9])+ theme_bw()
  p10 <- qplot(x=data[,10],geom="histogram") + geom_vline(xintercept=true[10],colour="red") + xlab(xlabels[10])+ theme_bw()
  p11 <- qplot(x=data[,11],geom="histogram") + geom_vline(xintercept=true[11],colour="red") + xlab(xlabels[11])+ theme_bw()
  p12 <- qplot(x=data[,12],geom="histogram") + geom_vline(xintercept=true[12],colour="red") + xlab(xlabels[12])+ theme_bw()
  p13 <- qplot(x=data[,13],geom="histogram") + geom_vline(xintercept=true[13],colour="red") + xlab(xlabels[13])+ theme_bw()
  p14 <- qplot(x=data[,14],geom="histogram") + geom_vline(xintercept=true[14],colour="red") + xlab(xlabels[14])+ theme_bw()
  p15 <- qplot(x=data[,15],geom="histogram") + geom_vline(xintercept=true[15],colour="red") + xlab(xlabels[15])+ theme_bw()
  p16 <- qplot(x=data[,16],geom="histogram") + geom_vline(xintercept=true[16],colour="red") + xlab(xlabels[16])+ theme_bw()
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16)
  #do.call(grid.arrange,p)
}

my_bs <- function (x, knots = NULL,  intercept = TRUE, degree = 3,
                   Boundary.knots = range(x)) 
{
  ord <- 1L + (degree <- as.integer(degree))
  if (ord <= 1) 
    stop("'degree' must be integer >= 1")
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax)) 
    x <- x[!nax]
  outside <- if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  }
  else FALSE
  
  Aknots <- sort(c(rep(Boundary.knots, ord), knots))
  if (any(outside)) {
    warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
    derivs <- 0:degree
    scalef <- gamma(1L:ord)
    basis <- array(0, c(length(x), length(Aknots) - degree - 
                          1L))
    e <- 1/4
    if (any(ol)) {
      k.pivot <- (1 - e) * Boundary.knots[1L] + e * Aknots[ord + 
                                                             1]
      xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree, 
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord, 
                         derivs)
      basis[ol, ] <- xl %*% (tt/scalef)
    }
    if (any(or)) {
      k.pivot <- (1 - e) * Boundary.knots[2L] + e * Aknots[length(Aknots) - 
                                                             ord]
      xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree, 
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord, 
                         derivs)
      basis[or, ] <- xr %*% (tt/scalef)
    }
    if (any(inside <- !outside)) 
      basis[inside, ] <- splineDesign(Aknots, x[inside], 
                                      ord)
  }
  else basis <- splineDesign(Aknots, x, ord)
  if (!intercept) 
    basis <- basis[, -1L, drop = FALSE]
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  #a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots, 
  #          Boundary.knots = Boundary.knots, intercept = intercept)
  #attributes(basis) <- c(attributes(basis), a)
  #class(basis) <- c("bs", "basis", "matrix")
  basis
}

my_qp <- function(X, dvec, amat, bvec,integer,true){
  return(solve.QP(X, dvec, amat, bvec,integer,true)$solution)
}

my_ns <- function (x, knots = NULL, intercept = FALSE, Boundary.knots = range(x)) 
{
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax)) 
    x <- x[!nax]
  if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
                                                   Boundary.knots[2L])
  }
  else outside <- FALSE

  nIknots <- length(knots)
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
  if (any(outside)) {
    basis <- array(0, c(length(x), nIknots + 4L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, x[ol] - k.pivot)
      tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0, 
                                                        1))
      basis[ol, ] <- xl %*% tt
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, x[or] - k.pivot)
      tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0, 
                                                        1))
      basis[or, ] <- xr %*% tt
    }
    if (any(inside <- !outside)) 
      basis[inside, ] <- splineDesign(Aknots, x[inside], 
                                      4)
  }
  else basis <- splineDesign(Aknots, x, 4)
  const <- splineDesign(Aknots, Boundary.knots, 4, c(2, 2))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L), 
                                                     drop = FALSE])
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  #a <- list(degree = 3L, knots = if (is.null(knots)) numeric() else knots, 
            #Boundary.knots = Boundary.knots, intercept = intercept)
  #attributes(basis) <- c(attributes(basis), a)
  #class(basis) <- c("ns", "basis", "matrix")
  basis
}

generate_data <- function(params,dist=1){
  #dist=1 Normal
  #    =2 skew normal
  #    =3 bimodal symmetric
  set.seed(1146999) #114, #1146, #11469, #1146999
  #number of samples
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
  mei1 <- 1500 + zig*zg + zib*zb + zia*za
  mei2 <- 2000 + zig*zg + zib*zb + zia*za
  mei3 <- 2300 + zig*zg + zib*zb + zia*za
  mei4 <- 2600 + zig*zg + zib*zb + zia*za
  mei5 <- 3200 + zig*zg + zib*zb + zia*za
  #sd EI
  sdei1 <- 40 +0#50
  sdei2 <- 40 +0#50
  sdei3 <- 80 +0#90
  sdei4 <- 40 +0#40
  sdei5 <- 60 +0#90
  
  #mean daily EE = mean daily EI
  mee1 <- 1500 + 130 + zeg*zg + zeb*zb + zea*za # 1500
  mee2 <- 2000 + 130 + zeg*zg + zeb*zb + zea*za # 2000
  mee3 <- 2200 + 130 + zeg*zg + zeb*zb + zea*za # 2300
  mee4 <- 2600 + 130 + zeg*zg + zeb*zb + zea*za # 2600
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
  df <- 8#5
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
  delta1 <- mvrnorm(n,rep(0,2),dmatrix)
  delta2 <- mvrnorm(n,rep(0,2),dmatrix)
  delta1es <- delta1[,1]# - delta1[,2]
  delta2es <- delta2[,1]# - delta2[,2]
  
  #calculate delta ES, positive change => ei > ee
  xes <- xei - xee
  
  #calculate true values T_ij of EE,EI, ES
  #   tei1 <- xei + delta1[,1]
  #   tei2 <- xei + delta2[,1]
  #   tee1 <- xee + delta1[,2]
  #   tee2 <- xee + delta2[,2]
  #   tes1 <- tei1 - tee1
  #   tes2 <- tei2 - tee2
  #calculate DLW est of EE
  if(dist==1){
    wee1 <- rnorm(n,xee,200) + delta1[,2] #truth 250 #rnorm(n,xee,.08*xee) + delta1[,2] #rnorm(n,xee,.05*rep(mee,each=n/2))
    wee2 <- rnorm(n,xee,200) + delta2[,2]#(n,xee,.08*xee) + delta2[,2]
    wee <- matrix(c(wee1,wee2),ncol=2,byrow=FALSE)
    #calculate DXA est of ES
    wes1 <- rnorm(n,xes,53) + delta1es#truth 72.862 #rnorm(n,xes,.04*abs(xes)+50) + delta1es
    wes2 <- rnorm(n,xes,53) + delta2es#rnorm(n,xes,.04*abs(xes)+50) + delta2es
    wes <- matrix(c(wes1,wes2),ncol=2,byrow=FALSE)
    #calc cheap EE
    yee1 <- be0+eecurve(xee)+geg*zg+geb*zb+gea*za + rnorm(n,0,380) + delta1[,2]#truth 408.534
    yee2 <- be0+eecurve(xee)+geg*zg+geb*zb+gea*za + rnorm(n,0,380) + delta2[,2]
    yee <- matrix(c(yee1,yee2),ncol=2,byrow=FALSE)
    #calc cheap ES
    yes1 <- bs0+escurve(xes)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) + delta1es#truth 216.1748
    yes2 <- bs0+escurve(xes)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) + delta2es
    yes <- matrix(c(yes1,yes2),ncol=2,byrow=FALSE)
  }
  if(dist==2){
    wee1 <- xee + rsnorm(n,0,200,10) + delta1[,2] #rnorm(n,xee,.05*rep(mee,each=n/2))
    wee2 <- xee + rsnorm(n,0,200,10) + delta2[,2]
    wee <- matrix(c(wee1,wee2),ncol=2,byrow=FALSE)
    #calculate DXA est of ES
    wes1 <- xes + rsnorm(n,0,53,10) + delta1es
    wes2 <- xes + rsnorm(n,0,53,10) + delta2es
    wes <- matrix(c(wes1,wes2),ncol=2,byrow=FALSE)
    #calc cheap EE
    yee1 <- be0+eecurve(xee)+geg*zg+geb*zb+gea*za + rsnorm(n,0,380,10) + delta1[,2]
    yee2 <- be0+eecurve(xee)+geg*zg+geb*zb+gea*za + rsnorm(n,0,380,10) + delta2[,2]
    yee <- matrix(c(yee1,yee2),ncol=2,byrow=FALSE)
    #calc cheap ES
    yes1 <- bs0+escurve(xes)+gig*zg+gib*zb+gia*za + rsnorm(n,0,210,10) + delta1es
    yes2 <- bs0+escurve(xes)+gig*zg+gib*zb+gia*za + rsnorm(n,0,210,10) + delta2es
    yes <- matrix(c(yes1,yes2),ncol=2,byrow=FALSE)
  }
  if(dist==3){
    wee1 <- xee + mix*rnorm(n,-175,sqrt(200^2-175^2)) + (1-mix)*rnorm(n,175,sqrt(200^2-175^2)) + delta1[,2]
    wee2 <- xee + mix*rnorm(n,-175,sqrt(200^2-175^2)) + (1-mix)*rnorm(n,175,sqrt(200^2-175^2)) + delta2[,2]
    wee <- matrix(c(wee1,wee2),ncol=2,byrow=FALSE)
    #calculate DXA est of ES
    wes1 <- xes + mix*rnorm(n,-45,sqrt(53^2-45^2)) + (1-mix)*rnorm(n,45,sqrt(53^2-45^2)) + delta1es
    wes2 <- xes + mix*rnorm(n,-45,sqrt(53^2-45^2)) + (1-mix)*rnorm(n,45,sqrt(53^2-45^2)) + delta2es
    wes <- matrix(c(wes1,wes2),ncol=2,byrow=FALSE)
    #calc cheap EE
    yee1 <- be0+eecurve(xee)+geg*zg+geb*zb+gea*za + mix*rnorm(n,-350,sqrt(380^2-350^2)) + (1-mix)*rnorm(n,350,sqrt(380^2-350^2)) + delta1[,2]
    yee2 <- be0+eecurve(xee)+geg*zg+geb*zb+gea*za + mix*rnorm(n,-350,sqrt(380^2-350^2)) + (1-mix)*rnorm(n,350,sqrt(380^2-350^2)) + delta2[,2]
    yee <- matrix(c(yee1,yee2),ncol=2,byrow=FALSE)
    #calc cheap ES
    yes1 <- bs0+escurve(xes)+gig*zg+gib*zb+gia*za + mix*rnorm(n,-190,sqrt(210^2-190^2)) + (1-mix)*rnorm(n,190,sqrt(210^2-190^2)) + delta1es
    yes2 <- bs0+escurve(xes)+gig*zg+gib*zb+gia*za + mix*rnorm(n,-190,sqrt(210^2-190^2)) + (1-mix)*rnorm(n,190,sqrt(210^2-190^2)) + delta2es
    yes <- matrix(c(yes1,yes2),ncol=2,byrow=FALSE)
  }
  return(list(xee=xee,xei=xei,xes=xes,wee=wee,wes=wes,yee=yee,yes=yes,zg=zg,zb=zb,za=za))
}

generate_data2 <- function(params,dist=1){
  #dist=1 Normal
  #    =2 skew normal
  #    =3 bimodal symmetric
  set.seed(1146999) #114, #1146, #11469, #1146999
  #number of samples
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
  mei1 <- 1500 + zig*zg + zib*zb + zia*za
  mei2 <- 2000 + zig*zg + zib*zb + zia*za
  mei3 <- 2300 + zig*zg + zib*zb + zia*za
  mei4 <- 2600 + zig*zg + zib*zb + zia*za
  mei5 <- 3200 + zig*zg + zib*zb + zia*za
  #sd EI
  sdei1 <- 40 +0#50
  sdei2 <- 40 +0#50
  sdei3 <- 80 +0#90
  sdei4 <- 40 +0#40
  sdei5 <- 60 +0#90
  
  #mean daily EE = mean daily EI
  mee1 <- 1500 + 130 + zeg*zg + zeb*zb + zea*za # 1500
  mee2 <- 2000 + 130 + zeg*zg + zeb*zb + zea*za # 2000
  mee3 <- 2200 + 130 + zeg*zg + zeb*zb + zea*za # 2300
  mee4 <- 2600 + 130 + zeg*zg + zeb*zb + zea*za # 2600
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
  df <- 8
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
  delta1 <- mvrnorm(n,rep(0,2),dmatrix)
  delta2 <- mvrnorm(n,rep(0,2),dmatrix)
  delta1es <- delta1[,1]# - delta1[,2]
  delta2es <- delta2[,1]# - delta2[,2]
  
  #calculate delta ES, positive change => ei > ee
  xes <- xei - xee
  
  #calculate true values T_ij of EE,EI, ES
  #   tei1 <- xei + delta1[,1]
  #   tei2 <- xei + delta2[,1]
  #   tee1 <- xee + delta1[,2]
  #   tee2 <- xee + delta2[,2]
  #   tes1 <- tei1 - tee1
  #   tes2 <- tei2 - tee2
  #calculate DLW est of EE
  if(dist==1){
    wee1 <- rnorm(n,xee+ delta1[,2],200)  #truth 250 #rnorm(n,xee,.08*xee) + delta1[,2] #rnorm(n,xee,.05*rep(mee,each=n/2))
    wee2 <- rnorm(n,xee+ delta2[,2],200) #(n,xee,.08*xee) + delta2[,2]
    wee <- matrix(c(wee1,wee2),ncol=2,byrow=FALSE)
    #calculate DXA est of ES
    wes1 <- rnorm(n,xes+ delta1es,53) #truth 72.862 #rnorm(n,xes,.04*abs(xes)+50) + delta1es
    wes2 <- rnorm(n,xes+ delta2es,53) #rnorm(n,xes,.04*abs(xes)+50) + delta2es
    wes <- matrix(c(wes1,wes2),ncol=2,byrow=FALSE)
    #calc cheap EE
    yee1 <- be0+eecurve(xee+ delta1[,2])+geg*zg+geb*zb+gea*za + rnorm(n,0,380) #truth 408.534
    yee2 <- be0+eecurve(xee+ delta2[,2])+geg*zg+geb*zb+gea*za + rnorm(n,0,380) 
    yee <- matrix(c(yee1,yee2),ncol=2,byrow=FALSE)
    #calc cheap ES
    yes1 <- bs0+escurve(xes+ delta1es)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) #truth 216.1748
    yes2 <- bs0+escurve(xes+ delta2es)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) 
    yes <- matrix(c(yes1,yes2),ncol=2,byrow=FALSE)
  }
  if(dist==2){
    wee1 <- xee + rsnorm(n,0,200,10) + delta1[,2] #rnorm(n,xee,.05*rep(mee,each=n/2))
    wee2 <- xee + rsnorm(n,0,200,10) + delta2[,2]
    wee <- matrix(c(wee1,wee2),ncol=2,byrow=FALSE)
    #calculate DXA est of ES
    wes1 <- xes + rsnorm(n,0,53,10) + delta1es
    wes2 <- xes + rsnorm(n,0,53,10) + delta2es
    wes <- matrix(c(wes1,wes2),ncol=2,byrow=FALSE)
    #calc cheap EE
    yee1 <- be0+eecurve(xee+ delta1[,2])+geg*zg+geb*zb+gea*za + rsnorm(n,0,380,10) 
    yee2 <- be0+eecurve(xee+ delta2[,2])+geg*zg+geb*zb+gea*za + rsnorm(n,0,380,10) 
    yee <- matrix(c(yee1,yee2),ncol=2,byrow=FALSE)
    #calc cheap ES
    yes1 <- bs0+escurve(xes+ delta1es)+gig*zg+gib*zb+gia*za + rsnorm(n,0,210,10) 
    yes2 <- bs0+escurve(xes+ delta2es)+gig*zg+gib*zb+gia*za + rsnorm(n,0,210,10) 
    yes <- matrix(c(yes1,yes2),ncol=2,byrow=FALSE)
  }
  if(dist==3){
    wee1 <- xee + mix*rnorm(n,-175,sqrt(200^2-175^2)) + (1-mix)*rnorm(n,175,sqrt(200^2-175^2)) + delta1[,2]
    wee2 <- xee + mix*rnorm(n,-175,sqrt(200^2-175^2)) + (1-mix)*rnorm(n,175,sqrt(200^2-175^2)) + delta2[,2]
    wee <- matrix(c(wee1,wee2),ncol=2,byrow=FALSE)
    #calculate DXA est of ES
    wes1 <- xes + mix*rnorm(n,-45,sqrt(53^2-45^2)) + (1-mix)*rnorm(n,45,sqrt(53^2-45^2)) + delta1es
    wes2 <- xes + mix*rnorm(n,-45,sqrt(53^2-45^2)) + (1-mix)*rnorm(n,45,sqrt(53^2-45^2)) + delta2es
    wes <- matrix(c(wes1,wes2),ncol=2,byrow=FALSE)
    #calc cheap EE
    yee1 <- be0+eecurve(xee+ delta1[,2])+geg*zg+geb*zb+gea*za + mix*rnorm(n,-350,sqrt(380^2-350^2)) + (1-mix)*rnorm(n,350,sqrt(380^2-350^2)) 
    yee2 <- be0+eecurve(xee + delta2[,2])+geg*zg+geb*zb+gea*za + mix*rnorm(n,-350,sqrt(380^2-350^2)) + (1-mix)*rnorm(n,350,sqrt(380^2-350^2))
    yee <- matrix(c(yee1,yee2),ncol=2,byrow=FALSE)
    #calc cheap ES
    yes1 <- bs0+escurve(xes+ delta1es)+gig*zg+gib*zb+gia*za + mix*rnorm(n,-190,sqrt(210^2-190^2)) + (1-mix)*rnorm(n,190,sqrt(210^2-190^2)) 
    yes2 <- bs0+escurve(xes+ delta2es)+gig*zg+gib*zb+gia*za + mix*rnorm(n,-190,sqrt(210^2-190^2)) + (1-mix)*rnorm(n,190,sqrt(210^2-190^2)) 
    yes <- matrix(c(yes1,yes2),ncol=2,byrow=FALSE)
  }
  return(list(xee=xee,xei=xei,xes=xes,wee=wee,wes=wes,yee=yee,yes=yes,zg=zg,zb=zb,za=za))
}

generate_data3 <- function(params,dist=1){
  #dist=1 Normal
  #    =2 skew normal
  #    =3 bimodal symmetric
  set.seed(1146999) #114, #1146, #11469, #1146999
  #number of samples
  n <- 500
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
  mei1 <- 1500 + zig*zg + zib*zb + zia*za
  mei2 <- 2100 + zig*zg + zib*zb + zia*za
  mei3 <- 2300 + zig*zg + zib*zb + zia*za
  mei4 <- 2900 + zig*zg + zib*zb + zia*za
  mei5 <- 3500 + zig*zg + zib*zb + zia*za
  #sd EI
  sdei1 <- 50 +0#50
  sdei2 <- 80 +0#50
  sdei3 <- 60 +0#90
  sdei4 <- 200 +0#40
  sdei5 <- 150 +0#90
  
  #mean daily EE = mean daily EI
  mee1 <- 1500 + 130 + zeg*zg + zeb*zb + zea*za # 1500
  mee2 <- 2100 + 130 + zeg*zg + zeb*zb + zea*za # 2000
  mee3 <- 2300 + 80 + zeg*zg + zeb*zb + zea*za # 2300
  mee4 <- 2900 + 130 + zeg*zg + zeb*zb + zea*za # 2600
  mee5 <- 3500 + 130 + zeg*zg + zeb*zb + zea*za #- 3500
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
  df <- 7
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
  delta1 <- mvrnorm(n,rep(0,2),dmatrix)
  delta2 <- mvrnorm(n,rep(0,2),dmatrix)
  delta3 <- mvrnorm(n,rep(0,2),dmatrix)
  delta4 <- mvrnorm(n,rep(0,2),dmatrix)
  delta5 <- mvrnorm(n,rep(0,2),dmatrix)
  delta6 <- mvrnorm(n,rep(0,2),dmatrix)
  delta1es <- delta1[,1]# - delta1[,2]
  delta2es <- delta2[,1]# - delta2[,2]
  delta3es <- delta3[,1]# - delta1[,2]
  delta4es <- delta4[,1]# - delta2[,2]
  delta5es <- delta5[,1]# - delta2[,2]
  delta6es <- delta6[,1]# - delta2[,2]
  
  
  #calculate delta ES, positive change => ei > ee
  xes <- xei - xee
  
  #calculate true values T_ij of EE,EI, ES
  #   tei1 <- xei + delta1[,1]
  #   tei2 <- xei + delta2[,1]
  #   tee1 <- xee + delta1[,2]
  #   tee2 <- xee + delta2[,2]
  #   tes1 <- tei1 - tee1
  #   tes2 <- tei2 - tee2
  #calculate DLW est of EE
  if(dist==1){
    wee1 <- rnorm(n,xee+ delta1[,2],200)  #truth 250 #rnorm(n,xee,.08*xee) + delta1[,2] #rnorm(n,xee,.05*rep(mee,each=n/2))
    wee2 <- rnorm(n,xee+ delta2[,2],200) #(n,xee,.08*xee) + delta2[,2]
    wee3 <- rnorm(n,xee+ delta3[,2],200)  #truth 250 #rnorm(n,xee,.08*xee) + delta1[,2] #rnorm(n,xee,.05*rep(mee,each=n/2))
    wee4 <- rnorm(n,xee+ delta4[,2],200) #(n,xee,.08*xee) + delta2[,2]
    wee5 <- rnorm(n,xee+ delta5[,2],200) #(n,xee,.08*xee) + delta2[,2]
    wee6 <- rnorm(n,xee+ delta6[,2],200) #(n,xee,.08*xee) + delta2[,2]
    
    wee <- matrix(c(wee1,wee2,wee3,wee4,wee5,wee6),ncol=6,byrow=FALSE)
    #calculate DXA est of ES
    wes1 <- rnorm(n,xes+ delta1es,53) #truth 72.862 #rnorm(n,xes,.04*abs(xes)+50) + delta1es
    wes2 <- rnorm(n,xes+ delta2es,53) #rnorm(n,xes,.04*abs(xes)+50) + delta2es
    wes3 <- rnorm(n,xes+ delta3es,53) #truth 72.862 #rnorm(n,xes,.04*abs(xes)+50) + delta1es
    wes4 <- rnorm(n,xes+ delta4es,53) #rnorm(n,xes,.04*abs(xes)+50) + delta2es
    wes5 <- rnorm(n,xes+ delta5es,53) #rnorm(n,xes,.04*abs(xes)+50) + delta2es
    wes6 <- rnorm(n,xes+ delta6es,53) #rnorm(n,xes,.04*abs(xes)+50) + delta2es
    
    wes <- matrix(c(wes1,wes2,wes3,wes4,wes5,wes6),ncol=6,byrow=FALSE)
    #calc cheap EE
    yee1 <- be0+eecurve(xee+ delta1[,2])+geg*zg+geb*zb+gea*za + rnorm(n,0,380) #truth 408.534
    yee2 <- be0+eecurve(xee+ delta2[,2])+geg*zg+geb*zb+gea*za + rnorm(n,0,380) 
    yee3 <- be0+eecurve(xee+ delta3[,2])+geg*zg+geb*zb+gea*za + rnorm(n,0,380) #truth 408.534
    yee4 <- be0+eecurve(xee+ delta4[,2])+geg*zg+geb*zb+gea*za + rnorm(n,0,380) 
    yee5 <- be0+eecurve(xee+ delta5[,2])+geg*zg+geb*zb+gea*za + rnorm(n,0,380) 
    yee6 <- be0+eecurve(xee+ delta6[,2])+geg*zg+geb*zb+gea*za + rnorm(n,0,380) 
    
    yee <- matrix(c(yee1,yee2,yee3,yee4,yee5,yee6),ncol=6,byrow=FALSE)
    #calc cheap ES
    yes1 <- bs0+escurve(xes+ delta1es)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) #truth 216.1748
    yes2 <- bs0+escurve(xes+ delta2es)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) 
    yes3 <- bs0+escurve(xes+ delta3es)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) #truth 216.1748
    yes4 <- bs0+escurve(xes+ delta4es)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) 
    yes5 <- bs0+escurve(xes+ delta5es)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) 
    yes6 <- bs0+escurve(xes+ delta6es)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) 
    
    yes <- matrix(c(yes1,yes2,yes3,yes4,yes5,yes6),ncol=6,byrow=FALSE)
  }
  if(dist==2){
    wee1 <- xee + rsnorm(n,0,200,10) + delta1[,2] #rnorm(n,xee,.05*rep(mee,each=n/2))
    wee2 <- xee + rsnorm(n,0,200,10) + delta2[,2]
    wee <- matrix(c(wee1,wee2),ncol=2,byrow=FALSE)
    #calculate DXA est of ES
    wes1 <- xes + rsnorm(n,0,53,10) + delta1es
    wes2 <- xes + rsnorm(n,0,53,10) + delta2es
    wes <- matrix(c(wes1,wes2),ncol=2,byrow=FALSE)
    #calc cheap EE
    yee1 <- be0+eecurve(xee+ delta1[,2])+geg*zg+geb*zb+gea*za + rsnorm(n,0,380,10) 
    yee2 <- be0+eecurve(xee+ delta2[,2])+geg*zg+geb*zb+gea*za + rsnorm(n,0,380,10) 
    yee <- matrix(c(yee1,yee2),ncol=2,byrow=FALSE)
    #calc cheap ES
    yes1 <- bs0+escurve(xes+ delta1es)+gig*zg+gib*zb+gia*za + rsnorm(n,0,210,10) 
    yes2 <- bs0+escurve(xes+ delta2es)+gig*zg+gib*zb+gia*za + rsnorm(n,0,210,10) 
    yes <- matrix(c(yes1,yes2),ncol=2,byrow=FALSE)
  }
  if(dist==3){
    wee1 <- xee + mix*rnorm(n,-175,sqrt(200^2-175^2)) + (1-mix)*rnorm(n,175,sqrt(200^2-175^2)) + delta1[,2]
    wee2 <- xee + mix*rnorm(n,-175,sqrt(200^2-175^2)) + (1-mix)*rnorm(n,175,sqrt(200^2-175^2)) + delta2[,2]
    wee <- matrix(c(wee1,wee2),ncol=2,byrow=FALSE)
    #calculate DXA est of ES
    wes1 <- xes + mix*rnorm(n,-45,sqrt(53^2-45^2)) + (1-mix)*rnorm(n,45,sqrt(53^2-45^2)) + delta1es
    wes2 <- xes + mix*rnorm(n,-45,sqrt(53^2-45^2)) + (1-mix)*rnorm(n,45,sqrt(53^2-45^2)) + delta2es
    wes <- matrix(c(wes1,wes2),ncol=2,byrow=FALSE)
    #calc cheap EE
    yee1 <- be0+eecurve(xee+ delta1[,2])+geg*zg+geb*zb+gea*za + mix*rnorm(n,-350,sqrt(380^2-350^2)) + (1-mix)*rnorm(n,350,sqrt(380^2-350^2)) 
    yee2 <- be0+eecurve(xee + delta2[,2])+geg*zg+geb*zb+gea*za + mix*rnorm(n,-350,sqrt(380^2-350^2)) + (1-mix)*rnorm(n,350,sqrt(380^2-350^2))
    yee <- matrix(c(yee1,yee2),ncol=2,byrow=FALSE)
    #calc cheap ES
    yes1 <- bs0+escurve(xes+ delta1es)+gig*zg+gib*zb+gia*za + mix*rnorm(n,-190,sqrt(210^2-190^2)) + (1-mix)*rnorm(n,190,sqrt(210^2-190^2)) 
    yes2 <- bs0+escurve(xes+ delta2es)+gig*zg+gib*zb+gia*za + mix*rnorm(n,-190,sqrt(210^2-190^2)) + (1-mix)*rnorm(n,190,sqrt(210^2-190^2)) 
    yes <- matrix(c(yes1,yes2),ncol=2,byrow=FALSE)
  }
  return(list(xee=xee,xei=xei,xes=xes,wee=wee,wes=wes,yee=yee,yes=yes,zg=zg,zb=zb,za=za))
}

generate_data4 <- function(params,dist=1,nrep=2){
  #dist=1 Normal
  #    =2 skew normal
  #    =3 bimodal symmetric
  set.seed(1146999) #114, #1146, #11469, #1146999
  #number of samples
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
  mei2 <- 2100 + zig*zg + zib*zb + zia*za
  mei3 <- 2300 + zig*zg + zib*zb + zia*za
  mei4 <- 2900 + zig*zg + zib*zb + zia*za
  mei5 <- 3200 + zig*zg + zib*zb + zia*za
  #sd EI
  sdei1 <- 50 +0#50
  sdei2 <- 80 +0#50
  sdei3 <- 60 +0#90
  sdei4 <- 150 +0#40
  sdei5 <- 220 +0#90
  
  #mean daily EE = mean daily EI
  mee1 <- 1500 + 130 + zeg*zg + zeb*zb + zea*za # 1500
  mee2 <- 2100 + 130 + zeg*zg + zeb*zb + zea*za # 2000
  mee3 <- 2300 + 80 + zeg*zg + zeb*zb + zea*za # 2300
  mee4 <- 2600 + 130 + zeg*zg + zeb*zb + zea*za # 2600
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
  delta <- list()
  for(i in 1:nrep){
    delta[[i]] = mvrnorm(n,rep(0,2),dmatrix)
  }

  
  #calculate delta ES, positive change => ei > ee
  #xes <- xei - xee
  
  xes <- sample(c(rnorm(n/5,-30,20),rnorm(n/5,30,20),rnorm(n/5,0,30),rnorm(n/5,-70,50),rnorm(n/5,70,50)),n,replace=FALSE)
  
  #calculate true values T_ij of EE,EI, ES
  #   tei1 <- xei + delta1[,1]
  #   tei2 <- xei + delta2[,1]
  #   tee1 <- xee + delta1[,2]
  #   tee2 <- xee + delta2[,2]
  #   tes1 <- tei1 - tee1
  #   tes2 <- tei2 - tee2
  #calculate DLW est of EE
  wee <- matrix(0,ncol=nrep,nrow=n,byrow=FALSE)
  wes <- matrix(0,ncol=nrep,nrow=n,byrow=FALSE)
  yee <- matrix(0,ncol=nrep,nrow=n,byrow=FALSE)
  yes <- matrix(0,ncol=nrep,nrow=n,byrow=FALSE)
  
  if(dist==1){
    for(i in 1:nrep){
      wee[,i] <- rnorm(n,xee+ delta[[i]][,2],200)  #truth 250 #rnorm(n,xee,.08*xee) + delta1[,2] #rnorm(n,xee,.05*rep(mee,each=n/2))
      #calculate DXA est of ES
      wes[,i] <- rnorm(n,xes+ delta[[i]][,1],53) #truth 72.862 #rnorm(n,xes,.04*abs(xes)+50) + delta1es
      #calc cheap EE
      yee[,i] <- be0+eecurve(xee+ delta[[i]][,2])+geg*zg+geb*zb+gea*za + rnorm(n,0,380) #truth 408.534
      #calc cheap ES
      yes[,i] <- bs0+escurve(xes+ delta[[i]][,1],k=0.14)+gig*zg+gib*zb+gia*za+ rnorm(n,0,210) #truth 216.1748
    }
  }
  if(dist==2){
    for(i in 1:nrep){
      wee[,i] <- xee + rsnorm(n,0,200,10) + delta[[i]][,2] #rnorm(n,xee,.05*rep(mee,each=n/2))
       #calculate DXA est of ES
      wes[,i] <- xes + rsnorm(n,0,53,10) + delta[[i]][,1]
      #calc cheap EE
      yee[,i] <- be0+eecurve(xee+ delta[[i]][,2])+geg*zg+geb*zb+gea*za + rsnorm(n,0,380,10) 
      #calc cheap ES
      yes[,i] <- bs0+escurve(xes+ delta[[i]][,1],k=0.1)+gig*zg+gib*zb+gia*za + rsnorm(n,0,210,10) 
    }
  }
  if(dist==3){
    for(i in 1:nrep){
      wee[,i] <- xee + mix*rnorm(n,-175,sqrt(200^2-175^2)) + (1-mix)*rnorm(n,175,sqrt(200^2-175^2)) + delta[[i]][,2]
       #calculate DXA est of ES
      wes[,i] <- xes + mix*rnorm(n,-45,sqrt(53^2-45^2)) + (1-mix)*rnorm(n,45,sqrt(53^2-45^2)) + delta[[i]][,1]
      #calc cheap EE
      yee[,i] <- be0+eecurve(xee+ delta[[i]][,2])+geg*zg+geb*zb+gea*za + mix*rnorm(n,-350,sqrt(380^2-350^2)) + (1-mix)*rnorm(n,350,sqrt(380^2-350^2)) 
      #calc cheap ES
      yes[,i] <- bs0+escurve(xes+ delta[[i]][,1],k=0.1)+gig*zg+gib*zb+gia*za + mix*rnorm(n,-190,sqrt(210^2-190^2)) + (1-mix)*rnorm(n,190,sqrt(210^2-190^2)) 
    }
  }
  return(list(xee=xee,xei=xei,xes=xes,wee=wee,wes=wes,yee=yee,yes=yes,zg=zg,zb=zb,za=za))
}
