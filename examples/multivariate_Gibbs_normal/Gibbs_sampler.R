#' Deterministic scanGibbs sampler for multivariate normal
#'
#' @param mu1 first component of mean parameter
#' @param mu2 second component of mean parameter
#' @param a variance of first component
#' @param b variance of second component
#' @param rho correlation between the two components
#' @param init initial value for Gibbs sampler
#' @param n size of Markov chain output
#'
#' @return Markov chain from DSGS for multivariate normal


Gibbs_sampler <- function(mu1, mu2, a, b, rho, init, n) {
    X <- matrix(0, nrow = n, ncol = 2)
    X[1, 1] = init[1]
    X[1, 2] =  init[2] 
    for (i in 2:n) {
        X[i, 1] = rnorm(1, mu1 + (rho / b) * (X[i - 1, 2] - mu2), sqrt(a - (rho ^ 2) / b))
        X[i, 2] = rnorm(1, mu2 + (rho / a) * (X[i, 1] - mu1), sqrt(b - (rho ^ 2) / a))
    }
    return(X)
}