rm(list=ls())
setwd("/Users/Bryan/Google Drive/Daniela/Microbiome")

source("getData.R")
source("XiaChenFungLiModel.R")

# n is how many top OTUs to keep
W <- W[,whichpart(colSums(W),n=75)]

set.seed(123)
#### Do EM
testOutput <- LNM.EM(W=W,base=50,EMiter=15,EMburn=5,MCiter=1500,MCburn=500,stepsize=0.01,p=0.05)

# take about 30 seconds each with above conditions
# save.image("microbiome826.RData")

## CHANGE: make LNM.EM output base
N <- nrow(W); Q <- ncol(W); base <- 50
M <- apply(W,1,sum)
eY <- tcrossprod(rep(1,N), testOutput$mu)
exp_eY <- exp(eY)
sum_exp_eY <- apply(exp_eY,1,sum)
eZ <- eW <- matrix(0,N,Q) 
for(i in 1:N){
  eZ[i,-base] <- exp_eY[i,]/(sum_exp_eY[i]+1)
  eZ[i,base] <- 1/(sum_exp_eY[i]+1)
  eW[i,] <- M[i]*eZ[i,]
}
# Mean squared prediction error
MSPE <- mean(apply((eW-W)^2,1,sum)/apply(W^2,1,sum)); round(MSPE,4)

  plot(as.vector(sqrt(W)),as.vector(sqrt(eW)),
       xlab=expression(sqrt(Observed)),
       ylab=expression(sqrt(Fitted)),main=paste("Fitting Graph"),
       pch=".",col="red")
  abline(a=0,b=1)
  mtext(paste( "(MSPE:",round(MSPE,3),")"))
  
# max one , sample 84, OTU 52
  W[84,52]
  # Why? Note that eY is the same for every OTU
  plot(W[,52])
  # Several negative outliers, will underfit
  inds <- which(W[,52]>1000); inds
  covars[inds,]
  table(covars$DayAmdmt)
  # Seems that it includes every DayAmdmt=12 points, and only 2 others
  
  plotdat <- data.frame(W=as.vector(sqrt(W)),eW=as.vector(sqrt(eW)),
                        DA=rep(covars$DayAmdmt,75),OTU=rep(colnames(W),each=119))
  
  plot(plotdat$W,plotdat$eW,xlab=expression(sqrt(Observed)),
       ylab=expression(sqrt(Fitted)),col=rainbow(9)[c(as.factor(plotdat$DA))],pch=20,
       main="Fitting Graph")
  legend("topright",sort(unique(as.character(plotdat$DA))),pch=20,col=rainbow(9)[sort(c(unique(as.factor(plotdat$DA))))])
  abline(a=0,b=1)
  mtext(paste( "(MSPE:",round(MSPE,3),")"))
  # library(lattice)
  # xyplot(plotdat$eW ~ plotdat$W, groups=plotdat$DA,
  #        key=list(text=list(as.character(unique(plotdat$DA))),
  #                 points=list(pch=20,col=as.factor(unique(plotdat$DA))),
  #                 corner=c(.95,1)),
  #        pch=20,col=1:8)
  # 
  # xyplot(plotdat$eW ~ plotdat$W, groups=plotdat$DA,
  #        auto.key=list(corner=c(.95,1)),
  #        pch=20)
  
  