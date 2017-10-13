#' Wsim
#'
#' Function to simulate raw counts W from model fit
#'
#' @param out model fit from LNM.EM
#' @param W raw count matrix used in model fit
#' @param niter number of simulations, defaults to 1000
#'
#' @export
Wsim <- function(out, W, niter = 1000) {
    N <- nrow(out$Y)
    Q <- ncol(out$Y) + 1
    base <- out$base
    W.m <- array(0, dim = c(N, Q, niter))
    M <- apply(W, 1, sum)
    for (i in 1:niter) {
        # set up Y take mu, sigma, simulate N Y_i's
        Y.m <- mvrnorm(n = N, mu = out$mu, Sigma = out$sigma)
        W.m[, , i] <- YtoW(Y = Y.m, M = M, base = base)
    }
    return(W.m)
}
