rm(list=ls())
setwd("/Users/Bryan/Google Drive/Daniela/Microbiome")

load("~/Google Drive/Daniela/Microbiome/output1121.RData")
plot(p1$Z,p1$eZ)
plot(p2$Z,p2$eZ)

getWindowVar <- function(out,window=0.01,slide=0.0005,expect=FALSE) {
  out <- data.frame(out)
  Z <- out[,1]
  eZ <- out[,2]

  endVal <- max(Z)
  # how many windows will we be looking at
  nWindow <- ceiling((endVal-window)/slide) + 1
  vals <- p_i <- c()
  for(i in 1:nWindow) {
    minVal <- 0 + (i-1)*slide
    maxVal <- window + (i-1)*slide
    subi <- Z[Z >= minVal & Z <= maxVal]
    if(expect==TRUE){
      subi <- eZ[Z >= minVal & Z <= maxVal]
    }
    if(length(subi)>1) {
      vals <- c(vals,var(subi))
      #p_i <- c(p_i,mean(c(minVal,maxVal)))
      p_i <- c(p_i,minVal)
    }
  }
  return(cbind(p_i,vals))
}


n11 <- sum(p1$Z)
n21 <- sum(p2$Z)
v11 <- getWindowVar(p1,window=0.05,slide=0.0051)
v21 <- getWindowVar(p2,slide=0.003)
plot(v11[,1],n11*v11[,1]*(1-v11[,1]))
points(v11)

pdf("movingWindowVarZ11.pdf")
par(mfrow=c(2,2))
#plot(getWindowVar(p1,window=0.1,slide=0.001),type="l",main="Window Size: 0.1",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p1,window=0.05,slide=0.001),type="l",main="Window Size: 0.05",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p1,window=0.025,slide=0.001),type="l",main="Window Size: 0.025",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p1,window=0.01,slide=0.001),type="l",main="Window Size: 0.01",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p1,window=0.005,slide=0.001),type="l",main="Window Size: 0.005",xlab=expression(p[i]),ylab="Within Window Variance")
par(mfrow=c(1,1))
dev.off();dev.off();

pdf("movingWindowVarZ21.pdf")
par(mfrow=c(2,2))
#plot(getWindowVar(p2,window=0.1,slide=0.001),type="l",main="Window Size: 0.1",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p2,window=0.05,slide=0.001),type="l",main="Window Size: 0.05",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p2,window=0.025,slide=0.001),type="l",main="Window Size: 0.025",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p2,window=0.01,slide=0.001),type="l",main="Window Size: 0.01",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p1,window=0.0005,slide=0.001),type="l",main="Window Size: 0.005",xlab=expression(p[i]),ylab="Within Window Variance")
par(mfrow=c(1,1))
dev.off();dev.off();

pdf("movingWindowVarEZ11.pdf")
par(mfrow=c(2,2))
#plot(getWindowVar(p1,window=0.1,slide=0.001,expect=TRUE),type="l",main="Window Size: 0.1",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p1,window=0.05,slide=0.001,expect=TRUE),type="l",main="Window Size: 0.05",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p1,window=0.025,slide=0.001,expect=TRUE),type="l",main="Window Size: 0.025",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p1,window=0.01,slide=0.001,expect=TRUE),type="l",main="Window Size: 0.01",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p1,window=0.0005,slide=0.001),type="l",main="Window Size: 0.005",xlab=expression(p[i]),ylab="Within Window Variance")
par(mfrow=c(1,1))
dev.off();dev.off();

pdf("movingWindowVarEZ21.pdf")
par(mfrow=c(2,2))
#plot(getWindowVar(p2,window=0.1,slide=0.001,expect=TRUE),type="l",main="Window Size: 0.1",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p2,window=0.05,slide=0.001,expect=TRUE),type="l",main="Window Size: 0.05",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p2,window=0.025,slide=0.001,expect=TRUE),type="l",main="Window Size: 0.025",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p2,window=0.01,slide=0.001,expect=TRUE),type="l",main="Window Size: 0.01",xlab=expression(p[i]),ylab="Within Window Variance")
plot(getWindowVar(p1,window=0.0005,slide=0.001),type="l",main="Window Size: 0.005",xlab=expression(p[i]),ylab="Within Window Variance")
par(mfrow=c(1,1))
dev.off();dev.off()




