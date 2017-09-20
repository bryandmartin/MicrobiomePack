rm(list=ls())
setwd("/Users/Bryan/Google Drive/Daniela/Microbiome")

source("getData.R")
source("XiaChenFungLiModel.R")

# Subset by experimental conditions
covars11 <- covars %>% filter(DayAmdmt == 11)
covars21 <- covars %>% filter(DayAmdmt == 21)
# get the samples of interest
samp11 <- covars11$X1
samp21 <- covars21$X1

# n is how many top OTUs to keep
W <- W[,whichpart(colSums(W),n=75)]

W11 <- W[which(S %in% samp11),]
W21 <- W[which(S %in% samp21),]


set.seed(123)
#### Do EM
out11 <- LNM.EM(W=W11,base=50,EMiter=10,EMburn=5,MCiter=1000,MCburn=500,
                stepsize=0.01,p=0.05,poorman=TRUE)
out21 <- LNM.EM(W=W21,base=50,EMiter=10,EMburn=5,MCiter=1000,MCburn=500,
                stepsize=0.01,p=0.05,poorman=TRUE)


compFitFun <- function(out,data) {
  N <- nrow(out$Y); Q <- ncol(out$Y)+1; base <- out$base
  eY <- tcrossprod(rep(1,N), out$mu)
  exp_eY <- exp(eY)
  sum_exp_eY <- apply(exp_eY,1,sum)
  eZ <- matrix(0,N,Q)
  for(i in 1:N){
    eZ[i,-base] <- exp_eY[i,]/(sum_exp_eY[i]+1)
    eZ[i,base] <- 1/(sum_exp_eY[i]+1)
  }
  
  Z <- makeComp(data)
  # Mean squared prediction error
  MSPE <- mean(apply((eZ-Z)^2,1,sum)/apply(Z^2,1,sum)); round(MSPE,4)
  
  plot(as.vector(Z),as.vector(eZ),
       xlab="Observed",
       ylab="Fitted",main=paste("Composition Fitting Graph"),
       pch=".",col="red")
  abline(a=0,b=1)
  mtext(paste( "(MSPE: ",round(MSPE,3),")",sep=""))
  return(list(Z=as.vector(Z),eZ=as.vector(eZ)))
}



p1 <- compFitFun(out11,W11); p2 <- compFitFun(out21,W21)
fullP <- cbind(p1$Z,p1$eZ)
fullP <- rbind(fullP, cbind(p2$Z,p2$eZ))
fullP <- cbind(fullP, rep(c(1,2),each=1200))
pdf("CombinedFittingGraph.pdf")
plot(fullP[,1],fullP[,2],
     xlab="Observed",
     ylab="Fitted",main=paste("Composition Fitting Graph"),
     pch=".",col=fullP[,3])
legend("bottomright",c("11","21"),col=c("black","red"),pch=c(19,19))
abline(a=0,b=1)
dev.off()


pdf("ExpectedCompDiff.pdf")
plot(unique(p1$eZ-p2$eZ),pch=20,ylab="Expected Difference",xlab="OTU Index",main="Expected Difference in Proportion by OTU")
dev.off()


# compare the variance of the mean estimates, with sigma
pdf("errorBars.pdf")
par(mfrow=c(2,1))
plot(out11$mu,pch=20,ylim=c(-5,5),ylab=expression(mu),main="DayAmdmt 11")
arrows(1:74, out11$mu-diag(out11$sigma), 1:74, out11$mu+diag(out11$sigma), length=0.05, angle=90, code=3)
arrows(1:74, out11$mu-apply(out11$mu.list,2,sd), 1:74, out11$mu+apply(out11$mu.list,2,sd), length=0.05, angle=90, code=3,col="red")
#arrows(1:74, out11$mu-apply(makeComp(W11)[,-50],2,sd), 1:74, out11$mu+apply(makeComp(W11)[,-50],2,sd), length=0.05, angle=90, code=3,col="red")


plot(out21$mu,pch=20,ylim=c(-5,5),ylab=expression(mu),main="DayAmdmt 21")
arrows(1:74, out21$mu-diag(out21$sigma), 1:74, out21$mu+diag(out21$sigma), length=0.05, angle=90, code=3)
arrows(1:74, out21$mu-apply(out21$mu.list,2,sd), 1:74, out21$mu+apply(out21$mu.list,2,sd), length=0.05, angle=90, code=3,col="red")
#arrows(1:74, out21$mu-apply(makeComp(W11)[,-50],2,sd), 1:74, out21$mu+apply(makeComp(W11)[,-50],2,sd), length=0.05, angle=90, code=3,col="red")

par(mfrow=c(1,1))
dev.off();dev.off()


