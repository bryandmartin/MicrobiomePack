rm(list = ls())
setwd("/Users/Bryan/Google Drive/Daniela/Microbiome")

source("getData.R")

# Subset by experimental conditions
covars11 <- covars %>% filter(DayAmdmt == 11)

# get the samples of interest
samp11 <- covars11$X1

# One trial: 50 most common + 250 random
# Another: 300 random
# n is how many top OTUs to keep
#m50 <- whichpart(colSums(W),n = 50)
#leftO <- seq(dim(W)[2])[-m50]
#r250 <- sample(leftO,300)
#W11 <- W[,c(r250)]
#W11 <- W[,sample(dim(W)[2],300)]
extractOTU <- function(x,nrand,nmost) {
  if (nmost >= 1) {
    mostI <- whichpart(colSums(x),n = nmost)
    leftO <- seq(dim(x)[2])[-mostI]
    randI <- sample(leftO,nrand)
    return(x[,c(mostI,randI)])
  } else {
    randI <- sample(seq(dim(x)[2]),nrand)
    return(x[,c(randI)])
  }
}
set.seed(123)
W11 <- extractOTU(W,300,0)

W11 <- W11[which(S %in% samp11),]

rownames(W11) <- as.character(1:nrow(W11))
require("DirichletMultinomial")
test <- list()
for (i in 1:4) {
  test[[i]] <- dmn(W11,i,verbose = TRUE)
}

lplc <- sapply(test,laplace)
plot(lplc,xlab = "Number of Dirichlet Components",ylab = "Model Fit")
best <- test[[which.min(lplc)]]; print(which.min(lplc))
# mixturewt(best)
# require("lattice")
# splom(log(fitted(best)))
# p0 <- fitted(test[[1]], scale = TRUE)     # scale by theta
# p4 <- fitted(best, scale = TRUE)
# colnames(p4) <- paste("m", 1:4, sep = "")
# (meandiff <- colSums(abs(p4 - as.vector(p0))))
# heatmapdmn(W11,test[[1]],test[[2]],30); dev.off()
# for now, only works with 1 mixture
# Just for fancy CM font plots, ignore
{
  require("extrafont")
  require("fontcm")
  # Takes a few
  # font_import()
  # font_install("fontcm")
  # loadfonts()
# Functions - TODO - wrap in file for sourcing
  savePDF <- function(filepath,img,latex=TRUE) {
    if (latex == TRUE) {
      pdf(filepath,family = "CM Roman")
      img
      dev.off() 
      embed_fonts(filepath,outfile = filepath)
    } 
    if (latex == FALSE) {
      pdf(filepath)
      img
      dev.off() 
    }
  }
}
dmmBars <- function(out,W,main="DMM Model Fit") {
  mu <- out@fit$Estimate / sum(out@fit$Estimate)
  eZ <- matrix(mu,byrow = TRUE,nrow = nrow(W), ncol = ncol(W))
  Z <- makeComp(W)
  par(mar = c(9,6,4,3) + .1,cex.main = 2)
  plot(as.vector(Z),as.vector(eZ),
       xlab = "",
       ylab = "",
       main = main,
       pch = 20)
  #abline(a=0,b=1)
  mtext(expression(frac(W[ij],sum(W[ij],j))),side = 1,line = 7.5,cex = 2)
  mtext("Scaled Dirichlet Component Estimates",side = 2,line = 3,cex = 1.5)
  eBarY <- eZ[1,]
  eBarX1 <-  out@fit$Lower / sum(out@fit$Lower)
  eBarX2 <-  out@fit$Upper / sum(out@fit$Upper)
  arrows(eBarX1,eBarY,eBarX2,eBarY,col = "red",code = 3,angle = 90,length = 0.01,lwd = 2)
}
dmmBarsT <- function(out,W,main="DMM Model Fit") {
  eZ <- makeComp(W)
  Q <- ncol(W); N <- nrow(W)
  
  Z <- rep(1:Q,each = N)
  par(mar = c(4.5,9,4,3) + .1,cex.main = 2)
  plot(as.vector(Z),as.vector(eZ),
       xlab = "",
       ylab = "",
       main = main,
       pch = 20)
  #abline(a=0,b=1)
  mtext("OTU Index",side = 1,line = 3,cex = 2)
  mtext(expression(frac(W[ij],sum(W[ij],j))),side = 2,las = 1,line = 3,cex = 2)
  eBarX <- 1:Q
  eBarY1 <-  out@fit$Lower / sum(out@fit$Lower)
  eBarY2 <-  out@fit$Upper / sum(out@fit$Upper)
  arrows(eBarX,eBarY1,eBarX,eBarY2,col = "red",code = 3,angle = 90,length = 0.01,lwd = 2)
}
savePDF("100817DMMPlots/DMM300.pdf",dmmBars(best,W11)); dev.off()
savePDF("100817DMMPlots/DMM300T.pdf",dmmBarsT(best,W11)); dev.off()
