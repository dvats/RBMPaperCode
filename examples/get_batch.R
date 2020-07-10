#' Divide data into batches
#'
#' @param X input data
#' @param n input dimensions
#' @param a number of batches
#' @param b batch size
#'
#' @return batched data
#' @export 


get_batch <- function(X, n, a, b) {
    p <- dim(X)[2]
    Z <- matrix(nrow = a, ncol = p)
    
    for (i in a:1) {
        Z[i,] = colMeans(X[(b * (i - 1) + 1):(b * i),])
    }

    return(Z)
}