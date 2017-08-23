######## FUNCTION ########
#### MC Step for a row
MCrow <- function(Yi,Wi,eYi,Q,base,sigma,MCiter,MCburn,stepsize=1) {
  Yi.MH <- matrix(0,MCiter,Q-1)
  aCount <- 0
  for(i in 1:MCiter) {
    # proposal
    Yi.star <- Yi + rnorm(Q-1,0,stepsize)
    # denominator
    Log_pt1 <- sum(Wi)*(log(sum(exp(Yi))+1)
                        -log(sum(exp(Yi.star))+1))
    # numerator
    Log_pt2 <- sum(Wi[-base]*(Yi.star-Yi))
    Log_pt3 <- -0.5*crossprod((Yi.star-eYi),solve(sigma))%*%(Yi.star-eYi)
    Log_pt4 <- -0.5*crossprod((Yi-eYi),solve(sigma))%*%(Yi-eYi)
    R_ratio <- Log_pt1+Log_pt2+Log_pt3-Log_pt4
    acceptance <- min(1,exp(R_ratio))
    temp <- runif(1)
    if(temp < acceptance) {
      Yi <- Yi.star
      aCount <- aCount + 1
    }
    Yi.MH[i,] <- Yi
  }
  Yi.new <- apply(Yi.MH[(MCburn+1):MCiter,],2,mean)
  ## returns a vector where first entry is acceptance rate, for ease of apply
  return(c(aCount/MCiter,Yi.new))
}

######## FUNCTION ########
#### Get MC matrix, acceptance
MCmat <- function(Y,W,eY,base,sigma,MCiter,MCburn,stepsize=1) {
  Q <- ncol(Y) + 1
  Y.MH <- sapply(1:nrow(W), function(i) MCrow(Yi=Y[i,],Wi=W[i,],eYi=eY[i,],base=base,sigma=sigma,MCiter=MCiter,MCburn=MCburn,stepsize=stepsize))
}


b0 <- attr(Y.p,"center")
View(rep(1,7769)%*%t(b0))
Center should be column means of purtubations
With no covariates, this is EY
# Equation 1, get additive logratio






save.image("microbiome.RData")
