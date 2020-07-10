#' average batch means
#'
#' @param Y Batched input data
#' @param m number of replications
#' @param n number of data points in each replication
#' @param a number of batches in each replication
#' @param b batch size
#'
#' @return Variance estimate using the average batch means method

#source("./../choosingBatch.R")
library(mcmcse)

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

ABM <- function(X, m, n, r, b) {
    ans1_r <- lapply(1:m, function(i)(mcse.multi(X[i,,], size = b, r = r, adjust = FALSE))$cov)
    abm_r <- (Reduce("+", ans1_r)) * (1 / m)

    # don't use lugsail if sum of eigen values < 0
    if(sum(diag(abm_r) < 0) > 0)
    {
        print("negative diag ABM")        
        abm <- ans1 <- lapply(1:m, function(i)(mcse.multi(X[i,,], size = b, r = 1, adjust = FALSE))$cov)
        return(abm)
    }

    # modify the covariance matrix if it is not positive definite
    min.eigen <- min(eigen(abm_r, only.values = TRUE)$values)
    if (min.eigen > 0) 
        return(abm_r)
    else {
        print("adjusting ABM")
        abm <- adjust_matrix(abm_r, n)
        return(abm)
    }
}