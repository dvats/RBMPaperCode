###########################################################
## Aim : Checking running estimates of ABM, naive method 
##       and variance estimator
##########################################################

set.seed(1)
source("./../get_batch.R")
source("./../ABM.R")
source("./../naive.R")
source("./../RBM.R")
source("./../ESS.R")
source("./../choosingbatch.R")

library(MCMCpack)
data(Nethvote)

running_multinomial <- function(rep, m, n, plot_step) {
    a <- floor(sqrt(n))
    b <- a
    # stores Markov chains
    X <- array(, dim = c(m, n, p))
    starting_buffer <- 500

    plot_mat_abm <- matrix(, nrow = rep, ncol = 1 + floor((n - starting_buffer) / plot_step))
    #plot_mat_nm <- matrix(, nrow = rep, ncol = floor((n - starting_buffer) / plot_step))
    plot_mat_rbm <- matrix(, nrow = rep, ncol = 1 + floor((n - starting_buffer) / plot_step))

    print(n + m)
    random <- matrix(sample(1:(rep * m), rep * m), nrow = rep)
    for (j in 1:rep) {
        print(j)
        #if ((j %% (rep / 10)) == 0) {
        #print(j)
        #}
        MLE <- MCMCmnl(vote ~ choicevar(distD66, "sqdist", "D66") +
                        choicevar(distPvdA, "sqdist", "PvdA") +
                        choicevar(distVVD, "sqdist", "VVD") +
                        choicevar(distCDA, "sqdist", "CDA") +
                        relig + class + income + educ + age + urban,
                        baseline = "D66", mcmc.method = "RWM", B0 = 0,
                        verbose = FALSE, mcmc = 1, thin = 1, tune = 0.5, burnin = 0, seed = j,
                        data = Nethvote)

        starting <- matrix(, nrow = m, ncol = p)
        for (i in 1:p) {
            starting[, i] <- seq(MLE[1, i] - 2, MLE[1, i] + 2, length = m)
        }

        # sample Markov chains
        for (i in 1:m) {

            X[i,,] <- MCMCmnl(vote ~ choicevar(distD66, "sqdist", "D66") +
            choicevar(distPvdA, "sqdist", "PvdA") +
            choicevar(distVVD, "sqdist", "VVD") +
            choicevar(distCDA, "sqdist", "CDA") +
            relig + class + income + educ + age + urban,
            baseline = "D66", beta.start = starting[i,], mcmc.method = "RWM", B0 = 0,
            verbose = 0, mcmc = n, thin = 1, tune = .5, burnin = 0, seed = i * j,
            data = Nethvote)
        }

        # calculate estimates across plot steps
        for (i in (0:floor((n - starting_buffer) / plot_step))) {
            b <- chooseBatch(X[, 1:(starting_buffer + plot_step * i),], r, 'sqroot')
            abm <- ABM(X[, 1:(starting_buffer + plot_step * i),], m, starting_buffer + plot_step * i, r, b)
            #nm <- NM(X[, (starting_buffer + 1):(starting_buffer + plot_step * i),], m, plot_step * i)
            rbm <- RBM(X[, 1:(starting_buffer + plot_step * i),], m, starting_buffer + plot_step * i, r, c, b)

            plot_mat_abm[j, i + 1] = (ESS(X[, 1:(starting_buffer + plot_step * i),], abm))[2]
            #plot_mat_nm[j, ((i - 1) * 5 + 1):(i * 5)] = c(nm[1, 1], nm[2, 2], nm[1, 2], det(nm), norm(nm, type = "F"))
            plot_mat_rbm[j, i + 1] = (ESS(X[, 1:(starting_buffer + plot_step * i),], rbm))[2]
        }

    }
    filename = paste("running_sqroot_2", "m", m, "n", n, "plot_step", plot_step, ".Rdata", sep = "_")
    save(plot_mat_abm, plot_mat_rbm, file = paste("output_files", filename, sep = "/"))

}


rep <- 1e3 
m_ <- c(2)
n <- 1e4
alpha <- 0.05
p <- 22
r <- 3 # Lugsail hyperparameters
c <- 0.5 # Lugsail hyperparameters
plot_step <- 100

for (iter1 in 1:length(m_)) {
    m <- m_[iter1]
    running_multinomial(rep, m, n, plot_step = plot_step)
}
