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