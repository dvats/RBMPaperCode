###########################################################
## Aim : Checking running estimates of ABM, naive method 
##       and variance estimator
##########################################################

set.seed(10)
library(mcmcse)
library(mcmc)
source("./../ABM.R")
source("./../naive.R")
source("./../RBM.R")
source("./../get_batch.R")
source("./../ESS.R")
source("./../choosingBatch.R")


# log unnormalized density
lufac <- function(A, B, mu) function(X) {
    r <- -A * (X[1] - mu) ^ 2 - B * (X[2] - X[1] ^ 2) ^ 2
    return(r)
}

running_rosenbrock <- function(rep, m, n, a, b, true_mean, A, B, mu) {
    print(n+m)
    a <- floor(sqrt(n))
    b <- a
    # stores Markov chains
    X <- array(, dim = c(m, n, 2))

    # log unnormalized posterior
    lupost = lufac(A, B, mu)

    coverage <- matrix(0, nrow = length(checkpoints) * rep, ncol = 3)
    determinant <- matrix(0, nrow = length(checkpoints) * rep, ncol = 3)
    ess_mat_concat <- matrix(, nrow = length(checkpoints) * rep, ncol = 3)
    ess_mat_avg <- matrix(, nrow = length(checkpoints) * rep, ncol = 3)

    ABM_estimate <- array(dim = c(length(checkpoints) * rep, 2, 2))
    NM_estimate <- array(dim = c(length(checkpoints) * rep, 2, 2))
    RBM_estimate <- array(dim = c(length(checkpoints) * rep, 2, 2))

    init <- cbind(seq(-edge1, edge1, length.out = m), sapply(seq(-edge2, edge2, length.out = m), function(t) rnorm(1, mean = t ^ 2, sd = sqrt(1 / (2 * B)))))

    for (j in 1:rep) {

        if ((j %% 100) == 0) 
            print(j)

        for (i in 1:m) {
            # MC is jumping rarely
            X[i,,] <- metrop(lupost, init[i,], n, scale = c(1, 3))$batch        ## Scale is 3 
        }

        # calculate covariance estimate at each checkpoint
        for (i in 1:length(checkpoints)) {
            b <- chooseBatch(X[, 1:checkpoints[i],], r, 'smart')
            ABM_estimate[(i - 1) * rep + j,,] <- ABM(X[, 1:checkpoints[i],], m, checkpoints[i], r, b)
            NM_estimate[(i - 1) * rep + j,,] <- NM(X[, 1:checkpoints[i],], m, checkpoints[i])
            RBM_estimate[(i - 1) * rep + j,,] <- RBM(X[, 1:checkpoints[i],], m, checkpoints[i], r, c, b)


            determinant[(i - 1) * rep + j, 1] = det(ABM_estimate[(i - 1) * rep + j,,])
            determinant[(i - 1) * rep + j, 2] = det(NM_estimate[(i - 1) * rep + j,,])
            determinant[(i - 1) * rep + j, 3] = det(RBM_estimate[(i - 1) * rep + j,,])

            # mean across all chains at given checkpoint
            M <- colMeans(X[, 1:checkpoints[i],], dims = 2)
            H <- qchisq(1 - alpha, p)

            # check if ABM is invertible
            is_invertible <- TRUE
            tryCatch({ result <- solve(ABM_estimate[(i - 1) * rep + j,,]) }, error = function(e) { is_invertible <<- FALSE })

            if (is_invertible) {
                if (m * checkpoints[i] * ((M - true_mean) %*% solve(ABM_estimate[(i - 1) * rep + j, , ]) %*% (M - true_mean)) < H) {
                    coverage[(i - 1) * rep + j, 1] = 1
                }
            }

            # check if NM is invertible
            is_invertible <- TRUE
            tryCatch({ result <- solve(NM_estimate[(i - 1) * rep + j,,]) }, error = function(e) { is_invertible <<- FALSE })

            if (is_invertible) {
                if (m * checkpoints[i] * ((M - true_mean) %*% solve(NM_estimate[(i - 1) * rep + j, , ]) %*% (M - true_mean)) < H) {
                    coverage[(i - 1) * rep + j, 2] = 1
                }
            }

            # check if RBM is invertible
            is_invertible <- TRUE
            tryCatch({ result <- solve(RBM_estimate[(i - 1) * rep + j,,]) }, error = function(e) { is_invertible <<- FALSE })

            if (is_invertible) {
                if (m * checkpoints[i] * ((M - true_mean) %*% solve(RBM_estimate[(i - 1) * rep + j, , ]) %*% (M - true_mean)) < H) {
                    coverage[(i - 1) * rep + j, 3] = 1
                }
            }


            ess_abm <- ESS(X[, 1:checkpoints[i],], ABM_estimate[(i - 1) * rep + j,,])
            ess_mat_concat[(i - 1) * rep + j, 1] <- ess_abm[1]
            ess_mat_avg[(i - 1) * rep + j, 1] <- ess_abm[2]

            ess_nm <- ESS(X[, 1:checkpoints[i],], NM_estimate[(i - 1) * rep + j,,])
            ess_mat_concat[(i - 1) * rep + j, 2] <- ess_nm[1]
            ess_mat_avg[(i - 1) * rep + j, 2] <- ess_nm[2]

            ess_rbm <- ESS(X[, 1:checkpoints[i],], RBM_estimate[(i - 1) * rep + j,,])
            ess_mat_concat[(i - 1) * rep + j, 3] <- ess_rbm[1]
            ess_mat_avg[(i - 1) * rep + j, 3] <- ess_rbm[2]

        }


    }
    filename = paste("CDE", edge1, edge2, "m", m, "n", n, ".Rdata", sep = "_")
    save(ABM_estimate, NM_estimate, RBM_estimate, ess_mat_concat, ess_mat_avg, determinant, coverage, file = paste("output_files", filename, sep = "/"))
}

edge1 = 8
edge2 = 8
rep <- 1e3 
m_ <- c(5, 10)
n <- 1e5
A = 1 / 20
B = 5
mu = 1
alpha <- 0.05
p <- 2
true_mean = c(mu, mu ^ 2 + 1 / (2 * A))
plot_step <- 500
checkpoints <- c(5e3, 1e4, 3e4, 5e4, 7e4, n)
r <- 3 # Lugsail hyperparameters
c <- 0.5 # Lugsail hyperparameters

for (i in 1:length(m_)) {
    m <- m_[i]
    running_rosenbrock(rep, m, n, a, b, true_mean, A, B, mu)
}


