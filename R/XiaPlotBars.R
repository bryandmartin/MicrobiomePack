#' XiaPlotBars
#'
#' Function to plot Xia model fit, with empirical 95% bars
#'
#' @param W raw count matrix used in model fit
#' @param out model fit from LNM.EM
#' @param main desired title for plot, defaults to Composition Fitting graph
#' @param niter number of simulations for bars, defaults to 1000
#' @param diagV boolean of whether to use diagonal Sigma for simulations, defaults to FALSE
#' @param smallV boolean of whether to use diagonal Sigma with very small entries for simulations, defaults to false
#' @param sub desired subtitle for plot, defaults to MSPE
#'
#' @export
XiaPlotBars <- function(W,out,main="Composition Fitting Graph",niter=1000,diagV=FALSE,smallV = FALSE,sub=FALSE) {
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
  mtext(expression(frac(W[ij],sum(W[ij],j))),side=1,line=7.5,cex=2)
  mtext(expression(phi^{-1}~(hat(mu)[j])),side=2,las=1,line=3,cex=2)
  
  if(diagV==TRUE) {
    temp <- matrix(0,Q-1,Q-1)
    diag(temp) <- diag(out$sigma)
    out$sigma <- temp
  }
  if(smallV==TRUE) {
    temp <- matrix(0,Q-1,Q-1)
    diag(temp) <- 1e-4
    out$sigma <- temp
  }
  if(sub!=FALSE) {
    mtext(sub)
  } else if (sub==FALSE) {
    mtext(paste( "(MSPE: ",round(MSPE,3),")",sep=""))
  }
  X.m <- Xsim(out,W,niter=niter)
  # place bars at every expected value (arbitrary choice 1)
  eBarY <- eZ[1,]
  eBarX1 <- apply(X.m,2,quantile,probs=c(0.025))
  eBarX2 <- apply(X.m,2,quantile,probs=c(0.975))
  arrows(eBarX1,eBarY,eBarX2,eBarY,col="red",code=3,angle=90,length=0.01)
}

#' XiaPlotBarsT
#'
#' Function to plot Xia model fit, with empirical 95% bars, transposed
#'
#' @param W raw count matrix used in model fit
#' @param out model fit from LNM.EM
#' @param main desired title for plot, defaults to Composition Fitting graph
#' @param niter number of simulations for bars, defaults to 1000
#' @param diagV boolean of whether to use diagonal Sigma for simulations, defaults to FALSE
#' @param smallV boolean of whether to use diagonal Sigma with very small entries for simulations, defaults to false
#' @param sub desired subtitle for plot, defaults to MSPE
#'
#' @export
XiaPlotBarsT <- function(W,out,main="Composition Fitting Graph",niter=1000,diagV=FALSE,smallV = FALSE,sub=FALSE) {
  N <- nrow(out$Y); Q <- ncol(out$Y)+1; base <- out$base
  eY <- tcrossprod(rep(1,N), out$mu)
  eZ <- YtoX(eY,base)
  Z <- makeComp(W)
  # Mean squared prediction error
  MSPE <- mean(apply((eZ-Z)^2,1,sum)/apply(Z^2,1,sum)); round(MSPE,4)
  par(mar=c(4.5,9,4,3)+.1,cex.main=2)
  eZ <- Z
  Z <- rep(1:Q,each=N)
  plot(as.vector(Z),as.vector(eZ),
       xlab="",
       ylab="",
       main=main,
       pch=20)
  #abline(a=0,b=1)
  mtext("OTU Index",side=1,line=3,cex=2)
  mtext(expression(frac(W[ij],sum(W[ij],j))),side=2,las=1,line=3,cex=2)
  if(diagV==TRUE) {
    temp <- matrix(0,Q-1,Q-1)
    diag(temp) <- diag(out$sigma)
    out$sigma <- temp
  }
  if(smallV==TRUE) {
    temp <- matrix(0,Q-1,Q-1)
    diag(temp) <- 1e-4
    out$sigma <- temp
  }
  if(sub!=FALSE) {
    mtext(sub)
  } else if (sub==FALSE) {
    mtext(paste( "(MSPE: ",round(MSPE,3),")",sep=""))
  }
  X.m <- Xsim(out,W,niter=niter)
  # place bars at every expected value (arbitrary choice 1)
  eBarX <- 1:Q
  eBarY1 <- apply(X.m,2,quantile,probs=c(0.025))
  eBarY2 <- apply(X.m,2,quantile,probs=c(0.975))
  arrows(eBarX,eBarY1,eBarX,eBarY2,col="red",code=3,angle=90,length=0.01)
}