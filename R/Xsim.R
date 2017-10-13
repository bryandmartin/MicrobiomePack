#' Xsim
#'
#' Function to simulate raw compositions X from model fit
#'
#' @param out model fit from LNM.EM
#' @param W raw count matrix used in model fit
#' @param niter number of simulations, defaults to 1000
#'
#' @export
Xsim <- function(out, W, niter = 1000) {
    N <- nrow(out$Y)
    Q <- ncol(out$Y) + 1
    base <- out$base
    X.m <- array(0, dim = c(N, Q, niter))
    M <- apply(W, 1, sum)
    for (i in 1:niter) {
        # set up Y take mu, sigma, simulate N Y_i's
        Y.m <- mvrnorm(n = N, mu = out$mu, Sigma = out$sigma)
        W.m <- YtoW(Y = Y.m, M = M, base = base)
        X.m[, , i] <- makeComp(W.m)
    }
    return(X.m)
}
