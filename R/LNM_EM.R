#' LNM.EM
#'
#' This function estimates the LNM model fit from Xia et al.
#' 
#' @param W count matrix, with OTUs as columns
#' @param X covariate matrix
#' @param base OTU index to be used for base
#' @param EMiter number of EM iterations, defaults to 10
#' @param EMburn number of EM iterations to burn, defaults to 5
#' @param MCiter number MC iterations, defaults to 1000
#' @param MCburn number of MC iterations to burn, defaults to 500
#' @param stepsize variance used for MH samples, defaults to 0.01. Tweak to adjust acceptance ratio
#' @param p size of purturbation used for logratios, defaults to 0.05
#' @param poorman boolean of whether to just use simple diagonal inverse to calculate sigma inverse. Defaults to FALSE
#'
#' @export
LNM.EM <- function(W, X, base, EMiter = 10, EMburn = 5, MCiter = 1000, MCburn = 500, stepsize = 0.01, 
    p = 0.05, poorman = FALSE) {
    # Base is value of D, stepsize is for MH, p is purturb Just take in W, calculate Y, eY, sigma Don't
    # need X yet, those are covariates
    W <- as.matrix(W)
    
    ### ROWS COLUMNS N is samples, Q is OTUs
    N <- nrow(W)
    Q <- ncol(W)
    X <- scale(X)
    # purturbed Y (N x Q-1) - function in getData.R
    Y.p <- logratios(W, base = base, p = p)
    # this is just column means (of OTUs)
    b0 <- attr(Y.p, "center")
    
    
    b <- OLS(X, Y.p)
    # This works only with no covariates Recall tcrossprod is x %*% t(y) N x 1 %*% 1 %*% Q eY is
    # constant across columns
    eY <- X %*% b + tcrossprod(rep(1, N), b0)
    
    # (Q-1) x (Q-1)
    sigma <- var(Y.p - eY)
    
    
    
    #### STORING RESULTS b0's by row of a matrix, rbind
    b0.list <- b0
    b.list <- b
    # sigma's by arraw, acomb3
    sigma.list <- sigma
    accept.list <- c()
    
    # Should be (MCiter x Q x N) (1000 x 75 x 119) Dont forget, first column is acceptance (ie want 74 x
    # 119 for data)
    for (em in 1:EMiter) {
        cat("EM iteration:", em, "\n")
        start <- proc.time()
        MCarray <- MCmat(Y = Y.p, W = W, eY = eY, N = N, Q = Q, base = base, sigma = sigma, MCiter = MCiter, 
            stepsize = stepsize, poorman = poorman)
        
        # should call 119 apply functions, each time, get 75 means. each of iteration values for OTU.
        # ORIGINAL for(i in 1:119) { Y.new[i,] <- apply(MCarray[(MCburn+1):MCiter,,i],2,mean) } ALTERNATIVE,
        # no for loop if large N, but transpose: seems faster by system.time
        Y.new <- t(apply(MCarray[(MCburn + 1):MCiter, , ], 3, colMeans))
        
        # recall first column
        accepts <- Y.new[, 1]
        Y.new <- Y.new[, 2:Q]
        
        # update beta, sigma Recall: b0 is means across OTUs
        b0 <- apply(Y.new, 2, mean)
        
        # sigSums <- matrix(0,Q-1,Q-1) # Get all sample matrices for(i in (MCburn+1):MCiter) { # first part
        # is Y sample Eps.samp <- t(MCarray[i,2:Q,1:N]) - eY sigSums <- sigSums + crossprod(Eps.samp) } CAN
        # EASILY MAKE ABOVE APPLY: next just take mean
        sigSumFun <- function(i) {
            return(crossprod(t(MCarray[i, 2:Q, 1:N]) - eY))
        }
        
        sigSum <- foreach(i = (MCburn + 1):MCiter, .combine = "+") %do% sigSumFun(i)
        
        sigma <- sigSum/(N * (MCiter - MCburn))
        
        b <- OLS(X, Y.new)
        
        eY <- X %*% b + tcrossprod(rep(1, N), b0)
        
        ### STORE after updating
        accept.list <- rbind(accept.list, accepts)
        b0.list <- rbind(b0.list, b0)
        b.list <- acomb3(b.list, b)
        sigma.list <- acomb3(sigma.list, sigma)
        end <- proc.time()
        cat(end[3] - start[3], "\n")
    }
    
    ## Next: take average of EM samples past burn. Included init, so have EMiter+1 total Note, below
    ## doesn't have any inital input
    accept.EM <- colMeans(accept.list[(EMburn):(EMiter), ])
    b0.EM <- colMeans(b0.list[(EMburn + 1):(EMiter + 1), ])
    b.EM <- apply(b.list[, , (EMburn + 1):(EMiter + 1)], c(1, 2), median)
    sigma.EM <- apply(sigma.list[, , (EMburn + 1):(EMiter + 1)], c(1, 2), mean)
    return(list(b0 = b0.EM, b = b.EM, sigma = sigma.EM, acceptance = accept.EM, b0.list = b0.list, b.list = b.list, 
        sigma.list = sigma.list, acceptance.list = accept.list, Y = Y.p, base = base))
}
