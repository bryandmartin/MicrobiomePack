require(doParallel)
require(abind)
require(doSNOW)
require(PDSCE)

######## FUNCTION ########
#### MC Step for a row
MCrow <- function(Yi,Wi,eYi,Q,base,sigInv,MCiter,stepsize=1) {
  # extra column for acceptance indicator
  Yi.MH <- matrix(0,MCiter,Q)
  for(i in 1:MCiter) {
    # proposal
    Yi.star <- Yi + rnorm(Q-1,0,stepsize)
    # denominator
    Eq5pt1 <- sum(Wi)*(log(sum(exp(Yi))+1)
                       -log(sum(exp(Yi.star))+1))
    # numerator
    Eq5pt2 <- sum(Wi[-base]*(Yi.star-Yi))
    Eq5pt3 <- -0.5*crossprod((Yi.star-eYi),sigInv)%*%(Yi.star-eYi)
    Eq5pt4 <- -0.5*crossprod((Yi-eYi),sigInv)%*%(Yi-eYi)
    fullRat <- Eq5pt1+Eq5pt2+Eq5pt3-Eq5pt4
    acceptance <- min(1,exp(fullRat))
    temp <- runif(1)
    aVal <- 0
    if(temp < acceptance) {
      Yi <- Yi.star
      aVal <- 1
    }
    # each row is a resampling of Yi, first col is whether accepted or not
    Yi.MH[i,] <- c(aVal,Yi)
  }
  
  ## returns a matrix
  #### Yi.MH: MH samples of row, first column ind if accepted (MCiter rows, Q cols)
  return(Yi.MH)
}


# helper
acomb3 <- function(...) abind(...,along=3)

######## FUNCTION ########
#### Get array
#### Dim 1 and 2 match MCrow
#### Dim 3 goes across rows of data (samples)
MCmat <- function(Y,W,eY,N,Q,base,sigma,MCiter,stepsize=1) {
  #sigInv <- solve(sigma)
  sigInv <- chol2inv(chol(sigma))
  
  MH_path <- function(i) {
    MCrow(Yi=Y[i,],Wi=W[i,],eYi=eY[i,],Q=Q,base=base,sigInv=sigInv,MCiter=MCiter,stepsize=stepsize)
  }
  #registerDoParallel(detectCores())
  Y.MH <-
    foreach(i=1:N,.combine='acomb3',.multicombine=TRUE) %do% {
      MH_path(i)
    }
  
  #stopImplicitCluster()
  # Should be (MCiter x Q x N)
  # Dont forget, first column is acceptance
  return(Y.MH)
}

#### Now have: MH samples of Y
#### Now need to: estimate parameters (M-step)
#### #### This requires: Calculations for the arguments of previous functions
#### #### And then: Plug in functions, follow formulas for estimates



######## FUNCTION ########
#### Do EM
LNM.EM <- function(W,base,EMiter=10,MCiter=1000,MCburn,stepsize=1,p=0.05) {
  # Base is value of D, stepsize is for MH, p is purturb
  # Just take in W, calculate Y, eY, sigma
  # Don't need X yet, those are covariates
  W <- as.matrix(W)
  
  ### ROWS COLUMNS
  # N is samples, Q is OTUs
  N <- nrow(W); Q <- ncol(W)
  # purturbed Y (N x Q-1) - function in getData.R
  Y.p <- logratios(W,base=base,p=p)
  # this is just column means (of OTUs)
  b0 <- attr(Y.p,"center")
  # This works only with no covariates
  # Recall tcrossprod is x %*% t(y)
  # N x 1 %*% 1 %*% Q
  # eY is constant across columns
  eY <- tcrossprod(rep(1,N), b0)
  
  # (Q-1) x (Q-1)
  sigma <- var(Y.p-eY)
  
  # Should be (MCiter x Q x N)
  # (1000 x 75 x 119)
  # Dont forget, first column is acceptance (ie want 74 x 119 for data)
  for(em in 1:EMiter) {
    MCarray <- MCmat(Y=Y.p,W=W,eY=eY,N=N,Q=Q,base=base,sigma=sigma,MCiter=MCiter,stepsize=stepsize)
    
    # should call 119 apply functions, each time, get 75 means. each of iteration values for OTU.
    # ORIGINAL
    # for(i in 1:119) {
    #   Y.new[i,] <- apply(MCarray[(MCburn+1):MCiter,,i],2,mean)
    # }
    # ALTERNATIVE, no for loop if large N, but transpose:
    # seems faster by system.time
    Y.new <- t(apply(MCarray[(MCburn+1):MCiter,,],3,colMeans))
    
    # recall first column
    accepts <- Y.new[,1]
    Y.new <- Y.new[,2:Q]
    
    # update beta, sigma
    # Recall: b0 is means across OTUs
    b0 <- apply(Y.new,2,mean)
    
    sigSums <- matrix(0,Q-1,Q-1)
    # Get all sample matrices
    for(i in (MHburn+1):MCiter) {
      # first part is Y sample
      Eps.samp <- t(MCarray[i,2:Q,1:N]) - eY
      sigSums <- sigSums + crossprod(Eps.samp)
    }
    ## CAN EASILY MAKE ABOVE APPLY: next just take mean
    
  }
  
  
}




# grab for debug:
base=50;EMiter=10;MCiter=100;MCburn=50;p=0.05;stepsize=0.01



# save.image("microbiome.RData")
