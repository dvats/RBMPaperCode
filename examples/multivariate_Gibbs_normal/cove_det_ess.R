###########################################################
## Aim : Checking running estimates of ABM, naive method 
##       and variance estimator
##########################################################

set.seed(10)
library(mcmcse)
source("Gibbs_sampler.R")
source("./../ABM.R")
source("./../naive.R")
source("./../RBM.R")
source("./../get_batch.R")
source("./../ESS.R")
source("./../choosingBatch.R")

running_mvg <- function(mu1, mu2, A, B, rho, true_mean, n, m, rep, plot_step) {
    print(n+m)
    print(rho)
    a <- floor(sqrt(n))
    b <- a
    # stores Markov chains
    X <- array(, dim = c(m, n, p))

    # specify method to calculate batch size
    if (rho == 0.5)
        method = 'sqroot'
    if (rho == 0.95)
        method = 'sqroot'
    if (rho == 0.999)
        method = 'smart'

    # true covariance matrix
    expected <- matrix(0, nrow = 2, ncol = 2)
    expected[1, 1] = A * ((A * B + rho ^ 2) / (A * B - rho ^ 2))
    expected[1, 2] = (2 * rho * A * B) / (A * B - rho ^ 2)
    expected[2, 1] = (2 * rho * A * B) / (A * B - rho ^ 2)
    expected[2, 2] = B * ((A * B + rho ^ 2) / (A * B - rho ^ 2))

    coverage <- matrix(0, nrow = length(checkpoints) * rep, ncol = 4)
    determinant <- matrix(0, nrow = length(checkpoints) * rep, ncol = 4)
    ess_mat_concat <- matrix(, nrow = length(checkpoints) * rep, ncol = 3)
    ess_mat_avg <- matrix(, nrow = length(checkpoints) * rep, ncol = 3)

    ABM_estimate <- array(dim = c(length(checkpoints) * rep, 2, 2))
    NM_estimate <- array(dim = c(length(checkpoints) * rep, 2, 2))
    RBM_estimate <- array(dim = c(length(checkpoints) * rep, 2, 2))

    # spread starting values
    starting <- matrix(0, nrow = m, ncol = 2)
    starting[, 1] <- seq(mu1 - 3 * A, mu1 + 3 * A, length = m)
    starting[, 2] <- seq(mu2 - 3 * B, mu2 + 3 * B, length = m)

    for (j in 1:rep) {

        if ((j %% 50) == 0) {
            print(j)
        }
         # sample Markov chains
        for (i in 1:m) {
            init <- starting[i,] # init is a vector now     #rnorm(1, mu1, sqrt(A))
            X[i,,] <- Gibbs_sampler(mu1, mu2, A, B, rho, init, n)
        }

        # calculate covariance estimate at each checkpoint
        for (i in 1:length(checkpoints)) {
            b <- chooseBatch(X[, 1:checkpoints[i],], r, method)
            ABM_estimate[(i - 1) * rep + j,,] <- ABM(X[, 1:checkpoints[i],], m, checkpoints[i], r, b)
            NM_estimate[(i - 1) * rep + j,,] <- NM(X[, 1:checkpoints[i],], m, checkpoints[i])
            RBM_estimate[(i - 1) * rep + j,,] <- RBM(X[, 1:checkpoints[i],], m, checkpoints[i], r, c, b)

        determinant[(i - 1) * rep + j, 1] = det(ABM_estimate[(i - 1) * rep + j,,])
        determinant[(i - 1) * rep + j, 2] = det(NM_estimate[(i - 1) * rep + j,,])
        determinant[(i - 1) * rep + j, 3] = det(RBM_estimate[(i - 1) * rep + j,,])
        determinant[(i - 1) * rep + j, 4] = det(expected)

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

        if (m * checkpoints[i] * ((M - c(mu1, mu2)) %*% solve(expected) %*% (M - c(mu1, mu2))) < H) {
            coverage[(i - 1) * rep + j, 4] = 1
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
    filename = paste("CDE", "rho", rho, "m", m, "n", n, ".Rdata", sep = "_")
    save(ABM_estimate, NM_estimate, RBM_estimate, ess_mat_concat, ess_mat_avg, determinant, coverage, file = paste("output_files", filename, sep = "/"))

}


mu1 <- 2
mu2 <- 50
A <- 1
B <- 1
n <- 1e4
rep <- 1e3
m_ <- c(10)
r <- 3
c <- 1 / 2
true_mean <- c(mu1, mu2)
plot_step <- 50
checkpoints <- c(1e2, 5e2, 1e3, 1.5e3, 2e3, 5e3, n)
p <- 2
alpha <- 0.05

for (i in 1:length(m_)) {
    m <- m_[i]

    rho = 0.5
    running_mvg(mu1, mu2, A, B, rho = .5, true_mean, n, m, rep, plot_step = plot_step)

    #rho = 0.95
    #running_mvg(mu1, mu2, A, B, rho = .95, true_mean, n, m, rep, plot_step = plot_step)

    rho = 0.999
    running_mvg(mu1, mu2, A, B, rho = .999, true_mean, n, m, rep, plot_step = plot_step)
}



