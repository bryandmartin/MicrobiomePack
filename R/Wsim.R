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
    mu <- get_mu(out)
    # if can sample from the same mu every time (no covariates)
    if (is.vector(mu)) {
        for (i in 1:niter) {
            # set up Y take mu, sigma, simulate N Y_i's
            Y.m <- mvrnorm(n = N, mu = mu, Sigma = out$sigma)
            W.m[, , i] <- YtoW(Y = Y.m, M = M, base = base)
        }
    }
    if (is.matrix(mu)) {
        for (i in 1:niter) {
            # apply out as vector stores as columns, transpose
            Y.m <- t(apply(mu, 1, function(x) mvrnorm(n = N, mu = x, Sigma = out$sigma)))
            W.m[, , i] <- YtoW(Y = Y.m, M = M, base = base)
        }
    }
    
    return(W.m)
}
