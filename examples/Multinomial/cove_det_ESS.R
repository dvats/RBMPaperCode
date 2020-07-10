########################################################################
## Aim : Calculate Variance using average batch means, 
##       replicated batch means and naive method for multinomial
########################################################################
source("./../get_batch.R")
source("./../ABM.R")
source("./../naive.R")
source("./../RBM.R")
source("./../ESS.R")
source("./../choosingbatch.R")

library(MCMCpack)
data(Nethvote)


running_multinomial <- function(rep, m, n, true_mean) {
    print(n + m)
    a <- floor(sqrt(n))
    b <- a
    # stores Markov chains
    X <- array(, dim = c(m, n, p))

    #determinant <- matrix(0, nrow = length(checkpoints) * rep, ncol = 3)
    ess_mat_concat <- matrix(, nrow = length(checkpoints) * rep, ncol = 3)
    ess_mat_avg <- matrix(, nrow = length(checkpoints) * rep, ncol = 3)

    ABM_estimate <- array(dim = c(length(checkpoints) * rep, p, p))
    NM_estimate <- array(dim = c(length(checkpoints) * rep, p, p))
    RBM_estimate <- array(dim = c(length(checkpoints) * rep, p, p))

    for (j in 1:rep) {

        if ((j %% (j/10)) == 0)
            print(j)

        # First draw of the Markov chain, 
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

        # calculate covariance estimate at each checkpoint
        for (i in 1:length(checkpoints)) {
            b <- chooseBatch(X[, 1:checkpoints[i],], r, 'sqroot')
            ABM_estimate[(i - 1) * rep + j,,] <- ABM(X[, 1:checkpoints[i],], m, checkpoints[i], r, b)
            NM_estimate[(i - 1) * rep + j,,] <- NM(X[, 1:checkpoints[i],], m, checkpoints[i])
            RBM_estimate[(i - 1) * rep + j,,] <- RBM(X[, 1:checkpoints[i],], m, checkpoints[i], r, c, b)


            #determinant[(i - 1) * rep + j, 1] = det(ABM_estimate[(i - 1) * rep + j,,])
            #determinant[(i - 1) * rep + j, 2] = det(NM_estimate[(i - 1) * rep + j,,])
            #determinant[(i - 1) * rep + j, 3] = det(RBM_estimate[(i - 1) * rep + j,,])

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
    filename = paste("CDE_sqroot_2", "m", m, "n", n, ".Rdata", sep = "_")
    save(ABM_estimate, NM_estimate, RBM_estimate, ess_mat_concat, ess_mat_avg, M, file = paste("output_files", filename, sep = "/"))

}

rep <- 1e3 
m<-2
m_ <- c(2)
n <- 1e4
alpha <- 0.05
p <- 22
checkpoints <- c(1e3, 3e3, 5e3, 7e3, n)
r <- 3 # Lugsail hyperparameters
c <- 0.5 # Lugsail hyperparameters


for (i in 1:length(m_)) {
    m <- m_[i]
    
    running_multinomial(rep, m, n, true_mean)
}
