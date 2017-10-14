#' XiaPlotWBars
#'
#' Function to plot Xia model fit for counts W, with empirical 95 percent points
#'
#' @param W raw count matrix used in model fit
#' @param out model fit from LNM.EM
#' @param main desired title for plot, defaults to Composition Fitting Graph
#' @param niter number of simulations for bars, defaults to 1000
#'
#' @export
XiaPlotWBars <- function(W, out, main = "Composition Fitting Graph", niter = 1000) {
    N <- nrow(out$Y)
    base <- out$base
    eY <- get_mu(out)
    if (is.vector(eY)) {
      eY <- tcrossprod(rep(1, N), eY)
    }
    M <- apply(W, 1, sum)
    eW <- YtoW(eY, M, base)
    # Mean squared prediction error
    MSPE <- mean(apply((eW - W)^2, 1, sum)/apply(W^2, 1, sum))
    round(MSPE, 4)
    par(mar = c(6, 12, 4, 6) + 0.1, cex.main = 2)
    plot(sqrt(as.vector(W)), sqrt(as.vector(eW)), xlab = "", ylab = "", main = main, pch = ".")
    # abline(a=0,b=1)
    mtext(paste("(MSPE: ", round(MSPE, 3), ")", sep = ""))
    mtext(expression(sqrt(W[ij])), side = 1, line = 4, cex = 2)
    mtext(expression(sqrt(n[i] ~ phi^{
        -1
    } ~ (hat(mu)[j]))), side = 2, las = 1, line = 2.5, cex = 2)
    
    W.m <- Wsim(out = out, W = W, niter = niter)
    # place bars at every expected value (arbitrary choice 1)
    ePointsL <- apply(W.m, c(1, 2), quantile, probs = c(0.025))
    ePointsH <- apply(W.m, c(1, 2), quantile, probs = c(0.975))
    points(sqrt(as.vector(W)), sqrt(as.vector(ePointsL)), col = "red", pch = ".")
    points(sqrt(as.vector(W)), sqrt(as.vector(ePointsH)), col = "red", pch = ".")
    points(sqrt(as.vector(W)), sqrt(as.vector(eW)), pch = ".")
    points(sqrt(as.vector(W))[as.vector(W) > as.vector(ePointsH)], sqrt(as.vector(eW))[as.vector(W) > 
        as.vector(ePointsH)], pch = 20)
    points(sqrt(as.vector(W))[as.vector(W) < as.vector(ePointsL)], sqrt(as.vector(eW))[as.vector(W) < 
        as.vector(ePointsL)], pch = 20)
}
