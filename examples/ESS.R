#' Calculate ESS
#'
#' @param X Markov chain output
#' @param cov_mat variance estimate 
#'
#' @return effective sample size

ESS <- function(X, cov_mat) {
    m <- dim(X)[1]
    n <- dim(X)[2]
    p <- dim(X)[3]
    Lambda = cov(X[1,,])
    concat_full <- X[1,,]
    for (i in 2:m) {
        concat_full <- rbind(concat_full, X[i,,])
        Lambda = Lambda + cov(X[i,,])
    }
    Lambda <- Lambda / m

    # ESS by concating the chains
    ess1 <- (det(cov(concat_full)) ^ (1 / p)) / (det(cov_mat) ^ (1 / p))
    det.Lambda <- eigen(Lambda, only.values = )

    # ESS by taking the avg of chains
    ess2 <-  exp ((determinant(Lambda, logarithm = TRUE)$modulus - determinant(cov_mat, logarithm = TRUE)$modulus)/p )
    return(c(ess1, ess2))
}