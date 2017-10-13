rm(list=ls())
setwd("/Users/Bryan/Google Drive/Daniela/Microbiome")

source("getData.R")
source("XiaChenFungLiModel.R")

# Subset by experimental conditions
covars11 <- covars %>% filter(DayAmdmt == 11)

# get the samples of interest
samp11 <- covars11$X1

# One trial: 50 most common + 250 random
# Another: 300 random
# n is how many top OTUs to keep
m50 <- whichpart(colSums(W),n=5)
leftO <- seq(dim(W)[2])[-m50]
r250 <- sample(leftO,295)
W11 <- W[,c(m50,r250)]
#W11 <- W[,sample(dim(W)[2],300)]
W11 <- W11[which(S %in% samp11),]



set.seed(123)
#### Do EM
out11 <- LNM.EM(W=W11,base=50,EMiter=10,EMburn=5,MCiter=1000,MCburn=500,
                stepsize=0.01,p=0.05,poorman=TRUE)

#savePDF("100517VarPlots/XiaPlotBars295v50.pdf",XiaPlotBars(W11,out11,main="Observed Data Fitting Graph",niter=niter,sub="(Full Sigma)"))

niter <- 1000
#savePDF("100517VarPlots/XiaPlotBars295v5.pdf",XiaPlotBars(W11,out11,main="Observed Data Fitting Graph",niter=niter,sub="(Full Sigma)"))
#savePDF("100517VarPlots/XiaPlotBars295v5T.pdf",XiaPlotBarsT(W11,out11,main="Observed Data Fitting Graph",niter=niter,sub="(Full Sigma)"))
#savePDF("100517VarPlots/XiaPlotBars300.pdf",XiaPlotBars(W11,out11,main="Observed Data Fitting Graph",niter=niter,sub="(Full Sigma)"))
#savePDF("100517VarPlots/XiaPlotBars300T.pdf",XiaPlotBarsT(W11,out11,main="Observed Data Fitting Graph",niter=niter,sub="(Full Sigma)"))


dev.off();dev.off()

