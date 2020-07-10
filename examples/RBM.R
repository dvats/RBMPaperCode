#' replicated batch means method
#'
#' @param Y Batched input data
#' @param m number of replications
#' @param n number of data points in each replication
#' @param a number of batches in each replication
#' @param b batch size
#'
#' @return Variance estimate using the replicated batch means method

source("./../get_batch.R")

# artificially inflate eigen values
adjust_matrix <- function(mat, N, epsilon = sqrt(log(N) / dim(mat)[2]), b = 9 / 10) {
    mat.adj <- mat
    adj <- epsilon * N ^ (-b)
    vars <- diag(mat)
    corr <- cov2cor(mat)
    eig <- eigen(corr)
    adj.eigs <- pmax(eig$values, adj)
    mat.adj <- diag(vars ^ (1 / 2)) %*% eig$vectors %*% diag(adj.eigs) %*% t(eig$vectors) %*% diag(vars ^ (1 / 2))
    return(mat.adj)
}

RBM <- function(X, m, n, r, c, b) {

    p <- dim(X)[3]
    if (b < 1 || b >= n || floor(n / b) <= 1)
        b <- floor(sqrt(n))
    a <- floor(n/b)
    b_r <- floor(b / r)
    a_r <- floor(n/b_r)
    Y1 <- matrix(0, nrow = m * a, ncol = p)
    Y1_r <- matrix(0, nrow = m * a_r, ncol = p)
    for (i in 1:m) {
        Y1[((i - 1) * a + 1):(i * a),] <- get_batch(X[i,,], n, a, b)
        Y1_r[((i - 1) * a_r + 1):(i * a_r),] <- get_batch(X[i,,], n, a_r, b_r)
    }
    overall_mean <- colMeans(Y1)
    ans3 <- lapply(1:(m * a), function(t) tcrossprod(Y1[t,] - overall_mean))
    ans3_r <- lapply(1:(m * a_r), function(t) tcrossprod(Y1_r[t,] - overall_mean))
    rbm <- (Reduce("+", ans3)) * (b / ((m * a) - 1))
    rbm_r <- (Reduce("+", ans3_r)) * (b_r / ((m * a_r) - 1))
    rbm_L <- (1 / (1 - c)) * rbm - (c / (1 - c)) * rbm_r

    # don't use lugsail if sum of eigen values < 0
    if(sum(diag(rbm_L) < 0) > 0)
    {
        print("negative diag RBM")         
        rtn <- rbm
        return(rtn)
    }

    # modify the covariance matrix if it is not positive definite
    min.eigen <- min(eigen(rbm_L, only.values = TRUE)$values)
    if(min.eigen  > 0)
        return(rbm_L)
    else{
        print("adjusting RBM")
        rbm_L = adjust_matrix(rbm_L, n)
        return(rbm_L)
    }

}