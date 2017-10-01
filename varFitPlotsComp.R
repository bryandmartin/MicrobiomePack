rm(list=ls())
setwd("/Users/Bryan/Google Drive/Daniela/Microbiome")

load("~/Google Drive/Daniela/Microbiome/output1121.RData")

# Just for fancy CM font plots, ignore
{
# require("extrafont")
# require("fontcm")
# Takes a few
# font_import()
# font_install("fontcm")
# loadfonts()
}

savePDF <- function(filepath,img,latex=TRUE) {
  if(latex==TRUE) {
    pdf(filepath,family="CM Roman")
    img
    dev.off() 
    embed_fonts(filepath,outfile=filepath)
  } 
  if(latex==FALSE) {
    pdf(filepath)
    img
    dev.off() 
  }
}

YtoX <- function(Y,base) {
  N <- nrow(Y); Q <- ncol(Y)+1
  exp_Y <- exp(Y)
  sum_exp_Y <- apply(exp_Y,1,sum)
  X <- matrix(0,N,Q)
  for(i in 1:N){
    X[i,-base] <- exp_Y[i,]/(sum_exp_Y[i]+1)
    X[i,base] <- 1/(sum_exp_Y[i]+1)
  }
  return(X)
}

YtoW <- function(Y,M,base) {
  N <- nrow(Y); Q <- ncol(Y)+1
  exp_Y <- exp(Y)
  sum_exp_Y <- apply(exp_Y,1,sum)
  X <- W <- matrix(0,N,Q)
  for(i in 1:N){
    X[i,-base] <- exp_Y[i,]/(sum_exp_Y[i]+1)
    X[i,base] <- 1/(sum_exp_Y[i]+1)
    W[i,] <- M[i]*X[i,]
  }
  return(W)
}

simpleMult <- function(W,main="Simple Multinomial Variance") {
  require("extrafont")
  W <- as.matrix(W)
  # W_ij / \sum_j W_{ij}
  # summing all column entries together? That's a rowSum
  Xaxis <- W/rowSums(W)
  # phat_j = \sum_i W_ij / \sum_ij W_ij
  phatj <- colSums(W) / sum(W)
  Yaxis <- matrix(phatj,nrow=nrow(W),ncol=ncol(W),byrow=TRUE)
  par(mar=c(9,9,4,3)+.1,cex.main=2)
  plot(Xaxis,Yaxis,pch=".",
       xlab="",
       ylab="",
       main=main)
  mtext(expression(frac(W[ij],sum(W[ij],j))),side=1,line=8,cex=2)
  mtext(expression(frac(sum(W[ij],i),sum(W[ij],ij))),side=2,las=1,line=3,cex=2)
  xx <- seq(0,max(makeComp(W))*1.1,length=100)
  lyy <- xx-2*sqrt(xx*(1-xx)/sum(W))
  uyy <- xx+2*sqrt(xx*(1-xx)/sum(W))
  lines(xx,lyy,col="red")
  lines(xx,uyy,col="red")
}

XiaPlot <- function(W,out,main="Composition Fitting Graph") {
  N <- nrow(out$Y); Q <- ncol(out$Y)+1; base <- out$base
  
  eY <- tcrossprod(rep(1,N), out$mu)
  eZ <- YtoX(eY,base)
  
  Z <- makeComp(W)
  # Mean squared prediction error
  MSPE <- mean(apply((eZ-Z)^2,1,sum)/apply(Z^2,1,sum)); round(MSPE,4)
  par(mar=c(6,6,4,1)+.1)
  plot(as.vector(Z),as.vector(eZ),
       xlab="",
       ylab="",
       main=main,
       pch=".")
  #abline(a=0,b=1)
  mtext(paste( "(MSPE: ",round(MSPE,3),")",sep=""))
  mtext(expression(frac(W[ij],sum(W[ij],j))),side=1,line=5)
  mtext(expression(phi^{-1}~(hat(mu)[j])),side=2,las=1,line=3)
}

XiaPlotW <- function(W,out,main="Observed Data Fitting Graph") {
  N <- nrow(out$Y); Q <- ncol(out$Y)+1; base <- out$base
  eY <- tcrossprod(rep(1,N), out$mu)
  M <- apply(W,1,sum)
  eW <- YtoW(eY,M,base)
  
  # Mean squared prediction error
  MSPE <- mean(apply((eW-W)^2,1,sum)/apply(W^2,1,sum)); round(MSPE,4)
  par(mar=c(5,7,4,1)+.1)
  plot(as.vector(W),as.vector(eW),
       xlab="",
       ylab="",
       main=main,
       pch=".")
  #abline(a=0,b=1)
  mtext(paste( "(MSPE: ",round(MSPE,3),")",sep=""))
  mtext(expression(W[ij]),side=1,line=3)
  mtext(expression(n[i]~phi^{-1}~(hat(mu)[j])),side=2,las=1,line=3)
}

# get simulations from model fit (W needed only to match counts)
Wsim <- function(out,W,niter=1000) {
  require(MASS)
  N <- nrow(out$Y); Q <- ncol(out$Y)+1; base <- out$base
  W.m <- array(0,dim=c(N,Q,niter))
  M <- apply(W,1,sum)
  for(i in 1:niter) {
    # set up Y
    # take mu, sigma, simulate N Y_i's
    Y.m <- mvrnorm(n=N,mu=out$mu,Sigma=out$sigma)
    W.m[,,i] <- YtoW(Y=Y.m,M=M,base=base)
  }
  return(W.m)
}

Xsim <- function(out,niter=1000) {
  require(MASS)
  N <- nrow(out$Y); Q <- ncol(out$Y)+1; base <- out$base
  X.m <- array(0,dim=c(N,Q,niter))
  for(i in 1:niter) {
    # set up Y
    # take mu, sigma, simulate N Y_i's
    Y.m <- mvrnorm(n=N,mu=out$mu,Sigma=out$sigma)
    X.m[,,i] <- YtoX(Y=Y.m,base=base)
  }
  return(X.m)
}


XiaPlotBars <- function(W,out,main="Composition Fitting Graph",niter=10000) {
  N <- nrow(out$Y); Q <- ncol(out$Y)+1; base <- out$base
  eY <- tcrossprod(rep(1,N), out$mu)
  eZ <- YtoX(eY,base)
  Z <- makeComp(W)
  # Mean squared prediction error
  MSPE <- mean(apply((eZ-Z)^2,1,sum)/apply(Z^2,1,sum)); round(MSPE,4)
  par(mar=c(9,10,4,3)+.1,cex.main=2)
  plot(as.vector(Z),as.vector(eZ),
       xlab="",
       ylab="",
       main=main,
       pch=20)
  #abline(a=0,b=1)
  mtext(paste( "(MSPE: ",round(MSPE,3),")",sep=""))
  mtext(expression(frac(W[ij],sum(W[ij],j))),side=1,line=7.5,cex=2)
  mtext(expression(phi^{-1}~(hat(mu)[j])),side=2,las=1,line=3,cex=2)
  
  X.m <- Xsim(out,niter=niter)
  # place bars at every expected value (arbitrary choice 1)
  eBarY <- eZ[1,]
  eBarX1 <- apply(X.m,2,quantile,probs=c(0.025))
  eBarX2 <- apply(X.m,2,quantile,probs=c(0.975))
  arrows(eBarX1,eBarY,eBarX2,eBarY,col="red",code=3,angle=90,length=0.01)
}

XiaPlotWBars <- function(W,out,main="Composition Fitting Graph",niter=10000) {
  N <- nrow(out$Y); Q <- ncol(out$Y)+1; base <- out$base
  eY <- tcrossprod(rep(1,N), out$mu)
  M <- apply(W,1,sum)
  eW <- YtoW(eY,M,base)
  # Mean squared prediction error
  MSPE <- mean(apply((eW-W)^2,1,sum)/apply(W^2,1,sum)); round(MSPE,4)
  par(mar=c(6,12,4,6)+.1,cex.main=2)
  plot(sqrt(as.vector(W)),sqrt(as.vector(eW)),
       xlab="",
       ylab="",
       main=main,
       pch=".")
  #abline(a=0,b=1)
  mtext(paste( "(MSPE: ",round(MSPE,3),")",sep=""))
  mtext(expression(sqrt(W[ij])),side=1,line=4,cex=2)
  mtext(expression(sqrt(n[i]~phi^{-1}~(hat(mu)[j]))),side=2,las=1,line=2.5,cex=2)
  
  W.m <- Wsim(out=out,W=W,niter=niter)
  # place bars at every expected value (arbitrary choice 1)
  ePointsL <- apply(W.m,c(1,2),quantile,probs=c(0.025))
  ePointsH <- apply(W.m,c(1,2),quantile,probs=c(0.975))
  points(sqrt(as.vector(W)),sqrt(as.vector(ePointsL)),col="red",pch=".")
  points(sqrt(as.vector(W)),sqrt(as.vector(ePointsH)),col="red",pch=".")
  points(sqrt(as.vector(W)),sqrt(as.vector(eW)),pch=".")
  points(sqrt(as.vector(W))[as.vector(W)>as.vector(ePointsH)],sqrt(as.vector(eW))[as.vector(W)>as.vector(ePointsH)],pch=20)
  points(sqrt(as.vector(W))[as.vector(W)<as.vector(ePointsL)],sqrt(as.vector(eW))[as.vector(W)<as.vector(ePointsL)],pch=20)
}


savePDF("093017VarPlots/simpleMult11.pdf",simpleMult(W11,main="Simple Multinomial Variance (11)"))
savePDF("093017VarPlots/simpleMult21.pdf",simpleMult(W21,main="Simple Multinomial Variance (21)"))

savePDF("093017VarPlots/XiaPlot11.pdf",XiaPlot(W11,out11,main="Composition Fitting Graph (11)"))
savePDF("093017VarPlots/XiaPlot21.pdf",XiaPlot(W21,out21,main="Composition Fitting Graph (21)"))


savePDF("093017VarPlots/XiaPlotW11.pdf",XiaPlotW(W11,out11,main="Observed Data Fitting Graph (11)"))
savePDF("093017VarPlots/XiaPlotW21.pdf",XiaPlotW(W21,out21,main="Observed Data Fitting Graph (21)"))

savePDF("093017VarPlots/XiaPlotBars11.pdf",XiaPlotBars(W11,out11,main="Composition Fitting Graph (11)"))
savePDF("093017VarPlots/XiaPlotBars21.pdf",XiaPlotBars(W21,out21,main="Composition Fitting Graph (21)"))


savePDF("093017VarPlots/XiaPlotWBars11.pdf",XiaPlotWBars(W11,out11,main="Observed Data Fitting Graph (11)"))
savePDF("093017VarPlots/XiaPlotWBars21.pdf",XiaPlotWBars(W21,out21,main="Observed Data Fitting Graph (21)"))

