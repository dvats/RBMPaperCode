###########################################################
## Aim : Checking running estimates of ABM, naive method 
##       and variance estimator
##########################################################

set.seed(666)
library(mcmcse)
source("Gibbs_sampler.R")
source("./../ABM.R")
source("./../naive.R")
source("./../RBM.R")
source("./../choosingBatch.R")

running_mvg <- function(mu1, mu2, A, B, rho, n, m, rep, plot_step)
{
    print(rho)
    # specify method to calculate batch size
    if (rho == 0.5) 
        method = 'sqroot'
    if (rho == 0.95)
        method = 'sqroot'
    if (rho == 0.999)
        method = 'smart'

    # stores Markov chains
    X <- array(, dim = c(m, n, 2))
    starting_buffer <- 0

    starting <- matrix(0, nrow = m, ncol = 2)
    starting[, 1] <- seq(mu1 - 3 * A, mu1 + 3 * A, length = m)
    starting[, 2] <- seq(mu2 - 3 * B, mu2 + 3 * B, length = m)

    plot_mat_abm <- matrix(, nrow = rep, ncol = 5 * floor((n - starting_buffer) / plot_step))
    plot_mat_nm <- matrix(, nrow = rep, ncol = 5 * floor((n - starting_buffer) / plot_step))
    plot_mat_rbm <- matrix(, nrow = rep, ncol = 5 * floor((n - starting_buffer) / plot_step))
    
    starting <- matrix(0, nrow = m, ncol = 2)
    starting[, 1] <- seq(mu1 - 3 * A, mu1 + 3 * A, length = m)
    starting[, 2] <- seq(mu2 - 3 * B, mu2 + 3 * B, length = m)
    print(n+m)
    for (j in 1:rep) {

        if ((j %% (rep/10)) == 0) {
            print(j)
        }

        # sample Markov chains
        for (i in 1:m) { 
            init <-  starting[i, ] 
            X[i,,] <- Gibbs_sampler(mu1, mu2, A, B, rho, init, n)
        }

        # calculate estimates across plot steps
        for (i in (1:floor((n - starting_buffer) / plot_step))) {
            b <- chooseBatch(X[, 1:(starting_buffer + plot_step * i),], r, method)

            abm <- ABM(X[, 1:(starting_buffer + plot_step * i),], m, starting_buffer + plot_step * i, r, b)
            nm <- NM(X[, 1:(starting_buffer + plot_step * i),], m, starting_buffer + plot_step * i)
            rbm <- RBM(X[, 1:(starting_buffer + plot_step * i),], m, starting_buffer + plot_step * i, r, c, b)

            plot_mat_abm[j, ((i - 1) * 5 + 1):(i * 5)] = c(abm[1, 1], abm[2, 2], abm[1, 2], det(abm), norm(abm, type = "F"))
            plot_mat_nm[j, ((i - 1) * 5 + 1):(i * 5)] = c(nm[1, 1], nm[2, 2], nm[1, 2], det(nm), norm(nm, type = "F"))
            plot_mat_rbm[j, ((i - 1) * 5 + 1):(i * 5)] = c(rbm[1, 1], rbm[2, 2], rbm[1, 2], det(rbm), norm(rbm, type = "F"))
        }

    }
    filename = paste("running", "rho", rho, "m", m, "n", n, ".Rdata", sep = "_")
    save(plot_mat_abm, plot_mat_nm, plot_mat_rbm, file = paste("output_files", filename, sep = "/"))
}


mu1 <- 2
mu2 <- 50
A <- 1
B <- 1
n <- 1e4
rep <- 1e2
m_ <- c(10)
r <- 3
c <- 1/2
plot_step <- 50
p <- 2
alpha <- 0.05



for (iter1 in 1:length(m_)) {
    m <- m_[iter1]
    running_mvg(mu1, mu2, A, B, rho = .5, n, m, rep, plot_step = plot_step)
    #running_mvg(mu1, mu2, A, B, rho = .95, n, m, rep, plot_step = plot_step)
    running_mvg(mu1, mu2, A, B, rho = .999, n, m, rep, plot_step = plot_step)

}


