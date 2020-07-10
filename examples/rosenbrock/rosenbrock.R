########################################################################
## Aim : Calculate Variance using average batch means, 
##       replicated batch means and stupid method for 2d Rosenbrock density
########################################################################
source("./../get_batch.R")
source("./../ABM.R")
source("./../stupid.R")
source("./../RBM.R")
source("./../choosingBatch.R")
source("./../ESS.R")

set.seed(11)
library(mcmc)
library(mcmcse)

# Log unnormalized posterior
lufac <- function(A, B, mu) function(X) {
    r <- -A * (X[1] - mu) ^ 2 - B * (X[2] - X[1] ^ 2) ^ 2
    return(r)
}

initial_stationary <-function(m, A, B, mu) {
    X <- matrix(, nrow = m, ncol = 2)
    X[, 1] <- rnorm(m, mean = mu, sd = sqrt(1 / (2 * A)))
    X[, 2] <- sapply(X[, 1], function(t) rnorm(n = 1, mean = t ^ 2, sd = sqrt(1 / (2 * B))))
    return(X)
}


running_rosenbrock <- function(rep, m, n, a, b, true_mean, A, B, mu, plot_step) {

    lupost = lufac(A, B, mu)
    starting_buffer = 0

    plot_mat_abm <- matrix(, nrow = rep, ncol = 1 + floor((n - starting_buffer) / plot_step))
    plot_mat_nm <- matrix(, nrow = rep, ncol = floor((n - starting_buffer) / plot_step))
    plot_mat_rbm <- matrix(, nrow = rep, ncol = 1 + floor((n - starting_buffer) / plot_step))

    # stores Markov chains
    X <- array(, dim = c(m, n, p))
    init <- cbind(seq(-edge1, edge1, length.out = m), sapply(seq(-edge2, edge2, length.out = m), function(t) rnorm(1, mean = t ^ 2, sd = sqrt(1 / (2 * B)))))

    print(n + m)
    for (j in 1:rep) {
        
        if ((j %% 10) == 0)
            print(j)

        # sample Markov chains
        for (i in 1:m) {
            # MC is jumping rarely
            X[i,,] <- metrop(lupost, initial = init[i,], nbatch = n, scale = c(1, 3))$batch
        }

        # calculate estimates across plot steps
        for (i in (1:floor((n - starting_buffer) / plot_step))) {
            if ((i %% 10) == 0)
                print(paste("Rep = ", j, ", % = ", i / (floor(n / plot_step))))
            b <- chooseBatch(X[, 1:(starting_buffer + plot_step * i),], r, 'smart')
            abm <- ABM(X[, 1:(starting_buffer + plot_step * i),], m, starting_buffer + plot_step * i, r, b)
            nm <- NM(X[, (starting_buffer + 1):(starting_buffer + plot_step * i),], m, plot_step * i)
            rbm <- RBM(X[, 1:(starting_buffer + plot_step * i),], m, starting_buffer + plot_step * i, r, c, b)

            plot_mat_abm[j, i] = (ESS(X[, 1:(starting_buffer + plot_step * i),], abm))[2]
            plot_mat_nm[j, i] = (ESS(X[, 1:(starting_buffer + plot_step * i),], nm))[2]
            plot_mat_rbm[j, i] = (ESS(X[, 1:(starting_buffer + plot_step * i),], rbm))[2]
        }

    }

    filename = paste("running_ESS", edge1, edge2, "m", m, "n", n, ".Rdata", sep = "_")
    save(plot_mat_abm, plot_mat_nm, plot_mat_rbm, file = paste("output_files", filename, sep = "/"))

}

edge1 = 8
edge2 = 8
rep <- 1e2 
m_ <- c(5, 10)
n_ <- c(1e5)
A = 1 / 20
B = 5
mu = 1
alpha <- 0.05
p <- 2
true_mean = c(mu, mu ^ 2 + 1 / (2 * A))
plot_step <- 500
r <- 3 # Lugsail hyperparameters
c <- 0.5 # Lugsail hyperparameters

for (iter1 in 1:length(m_)) {
    for (iter2 in 1:length(n_)) {
        m <- m_[iter1]
        n <- n_[iter2]
        a <- floor(sqrt(n))
        b <- a
        paste0("m is ",m, " n is ", n)
        running_rosenbrock(rep, m, n, a, b, true_mean, A, B, mu, plot_step)

    }
}
