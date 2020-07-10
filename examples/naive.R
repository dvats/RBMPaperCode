#' naive method
#'
#' @param Y Batched input data
#' @param m number of replications
#' @param n number of data points in each replication
#' @param a number of batches in each replication
#' @param b batch size
#'
#' @return Variance estimate using the naive method

NM <- function(X, m, n) {
    a <- floor(sqrt(n))
    b <- a
    p <- dim(X)[3]
    Ymean <- matrix(, nrow = m, ncol = p)
    for (i in 1:m) {
        Ymean[i,] = colMeans(X[i,,])
    }
    nm <- n * var(Ymean)
    return(nm)
}