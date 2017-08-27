rm(list=ls())
setwd("/Users/Bryan/Google Drive/Daniela/Microbiome")

source("getData.R")

# n is how many top OTUs to keep
W <- W[,whichpart(colSums(W),n=75)]

# Assume flat dirichlet (marginally uniform)
alpha.prior <- rep(1,75)

# Multinomial, draw N things, put them into K boxes based on probability
alpha.post <- alpha.prior + colSums(W)

# mean vector for probability
mu.post <- alpha.post/sum(alpha.post)*sum(W)

# covariance for probability
getDMvar <- function(n,ai,asum) {
  return(n * (ai/asum) * (1-ai/asum) * ((n+asum)/(1+asum)))
}
getDMcovar <- function(n,ai,aj,asum) {
  return(-n * ((ai*aj)/(asum^2)) * ((n+asum)/(1+asum)) )
}
K <- ncol(W)
sig.post <- matrix(0,nrow=K,ncol=K)
for(i in 1:K) {
  for(j in 1:K) {
    if(i==j) {
      sig.post[i,j] <- getDMvar(n=sum(W),ai=alpha.post[i],asum=sum(alpha.post))
    } else {
      sig.post[i,j] <- getDMcovar(n=sum(W),ai=alpha.post[i],aj=alpha.post[j],asum=sum(alpha.post))
    }
  }
}

library(lattice)
diag(sig.post) <- 0
levelplot(sig.post)
levelplot(testOutput$sigma)
# check dirmult package? different results than:
# temp <- dirmult(W)
# Even when using init = rep(1,75)
# Why is it converging?