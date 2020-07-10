#' average batch means
#'
#' @param X array of m by n by p chains
#' @param r lugsail parameter
#' @return average batch size across all chains

library(mcmcse)
chooseBatch <- function(X, r, method = 'sqroot') 
{
	p <- dim(X)[3]
    n <- dim(X)[2]
    if (method == 'sqroot') 
        b <- floor(sqrt(n))
    if (method == 'smart')
        b <- floor(mean(sapply(1:m, function(t) batchSize(X[t,,]))))
    b <- min(b, floor(n / (p + 1))) # making sure covariance matrix won't be necessarily non positive definite
    b <- max(b, 2 * r) # making sure lugsail works
    return(b)
}