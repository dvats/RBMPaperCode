######################################
##### MVG
######################################
set.seed(1)
library(gridExtra)
library(scales)
library(plotly)
library(mcmc)
library(MCMCpack)
library(plot3D)

## TRACE PLOTS
# Sample Markov chain from Gibbs sampler
Gibbs_sampler <- function(mu1, mu2, a, b, rho, init, n) {
    X <- matrix(0, nrow = n, ncol = 2)
    X[1, 1] = init[1]
    X[1, 2] = init[2] 
    for (i in 2:n) {
        X[i, 1] = rnorm(1, mu1 + (rho / b) * (X[i - 1, 2] - mu2), sqrt(a - (rho ^ 2) / b))
        X[i, 2] = rnorm(1, mu2 + (rho / a) * (X[i, 1] - mu1), sqrt(b - (rho ^ 2) / a))
    }
    return(X)
}

make_trace_plots <- function(m, n, p, A, B, rho) {

    starting <- matrix(0, nrow = m, ncol = 2) 
    starting[, 1] <- seq(mu1 - 3 * A, mu1 + 3 * A, length = m)
    starting[, 2] <- seq(mu2 - 3 * B, mu2 + 3 * B, length = m)

    X <- array(, dim = c(m, n, p)) # stores the Markov chain
    for (i in 1:m) {
        init <- starting[i,]
        X[i,,] <- Gibbs_sampler(mu1, mu2, A, B, rho, init, n)
    }

    plot_name <- paste("MVG_trace_plots", rho, ".pdf", sep = "_")
    pdf(plot_name, height = 6, width = 6)
    yrange = range(X[1,, 1], X[2,, 1])
    plot.ts(X[1,, 1], col = alpha("blue", 1), type = 'l', lty = 1, ylim = c(yrange[1] / 1.05, yrange[2] * 1.05), ylab = 'component 1', xlab = 'Chain Length', cex.axis = 1.4, cex.lab = 1.4)
    lines(X[2,, 1], col = alpha("red", 0.6), lty = 2)

    op <- par(cex = 1.2)
    legend("topright", legend = c("chain 1", "chain 2"), col = c("blue", "red"), lty = c(1, 2), bty = 'n')

    dev.off()

}

mu1 <- 2
mu2 <- 50
A <- 1
B <- 1
n <- 1e4
m <- 2
p <- 2
rho <- 0.5
make_trace_plots(m, n, p, A, B, rho)

rho <- 0.999
make_trace_plots(m, n, p, A, B, rho)


## PLOT FUNCTION FOR RUNNING PLOTS
make_plots <- function(A, B, rho, m, n, rep, plot_step = 50) {

    filename = paste("MVG_running", "rho", rho, "m", m, "n", n, ".Rdata", sep = "_")
    load(filename)
    starting_buffer <- 0
    # expected is the true value of covariance matrix
    expected <- matrix(0, nrow = 2, ncol = 2)
    expected[1, 1] = A * ((A * B + rho ^ 2) / (A * B - rho ^ 2))
    expected[1, 2] = (2 * rho * A * B) / (A * B - rho ^ 2)
    expected[2, 1] = (2 * rho * A * B) / (A * B - rho ^ 2)
    expected[2, 2] = B * ((A * B + rho ^ 2) / (A * B - rho ^ 2))

    if (rho == 0.5)
        legend_loc = "topright"
    if (rho == 0.999)
        legend_loc = "bottomright"

    truths <- c(expected[1, 1], expected[2, 2], expected[1, 2], det(expected), norm(expected, "F"))
    names <- c("Var1", "Var2", "Covariance", "Determinant", "Frobenius Norm")
    for (k in 5:5) {
        ind <- (0:(floor((n - starting_buffer) / plot_step) - 1)) * 5 + k
        x <- (0:(floor((n - starting_buffer) / plot_step) - 1)) * plot_step + plot_step
        plot_name <- paste("MVG_running", "rho", rho, "m", m, "n", n, "type", k, ".pdf", sep = "_")
        pdf(paste(plot_name, sep = "/"), height = 6, width = 6)

        yrange = range(colMeans(plot_mat_nm[, ind]), colMeans(plot_mat_abm[, ind]), colMeans(plot_mat_rbm[, ind]))
        plot(x, colMeans(plot_mat_abm[, ind]), type = "l", col = "green", lty = 1, lwd = 2, ylim = c(yrange[1]/1.15, 1.05*yrange[2]), ylab = names[k], xlab = "Chain Length", cex.lab = 1.4, cex.axis = 1.4)
        se_abm = sqrt(apply(plot_mat_abm[, ind], 2, var) / dim(plot_mat_abm)[1])
        segments(x0 = x, y0 = colMeans(plot_mat_abm[, ind]) - 1.96 * se_abm, x1 = x, y1 = colMeans(plot_mat_abm[, ind]) + 1.96 * se_abm, col = alpha("green", 0.8))

        lines(x, colMeans(plot_mat_nm[, ind]), type = "l", col = "red", lty = 2, lwd = 2)
        se_sm = sqrt(apply(plot_mat_nm[, ind], 2, var) / dim(plot_mat_nm)[1])
        segments(x0 = x, y0 = colMeans(plot_mat_nm[, ind]) - 1.96 * se_sm, x1 = x, y1 = colMeans(plot_mat_nm[, ind]) + 1.96 * se_sm, col = alpha("red", 0.3))

        lines(x, colMeans(plot_mat_rbm[, ind]), type = "l", col = "blue", lty = 6, lwd = 2)
        se_rbm = sqrt(apply(plot_mat_rbm[, ind], 2, var) / dim(plot_mat_rbm)[1])
        segments(x0 = x, y0 = colMeans(plot_mat_rbm[, ind]) - 1.96 * se_rbm, x1 = x, y1 = colMeans(plot_mat_rbm[, ind]) + 1.96 * se_rbm, col = alpha("blue", 0.4))

        abline(h = truths[k], lty = 2)
        op <- par(cex = 1.2)
        legend(legend_loc, legend = c("ABM", "Naive", "RBM"), col = c("green", "red", "blue"), lty = c(1, 2, 6), bty = 'n')
        dev.off()
    }
}

A <- 1
B <- 1
n <- 1e4
m_ <- c(5,10)
rep <- 1e2
plot_step <- 50

for (m in m_) {
    rho <- 0.5
    make_plots(A, B, rho, m, n, rep, plot_step = plot_step)

    rho <- 0.999
    make_plots(A, B, rho, m, n, rep, plot_step = plot_step)
}


## MAKE TABLES
pluck <- function(mat, k) {
    if (k == 1) {
        foo <- mat[, 1, 1]
        return(foo)
    }
    if (k == 2) {
        foo <- mat[, 2, 2]
        return(foo)
    }
    if (k == 3) {
        foo <- mat[, 1, 2]
        return(foo)
    }
    if (k == 4) {
        foo <- apply(mat, 1, det)
        return(foo)
    }

    if (k == 5) {
        foo <- apply(mat, 1, norm, "F")
        return(foo)
    }
}

make_plots <- function(A, B, rho, m, n, rep, plot_step = 50, type) {

    filename = paste("MVG_CDE", "rho", rho, "m", m, "n", n, ".Rdata", sep = "_")
    load(filename)

    # Tables:
    if (type == 'coverage') {

        coverage = rbind(coverage[1:(rep * 3),], coverage[(6 * rep + 1):(7 * rep),])
        determinant = rbind(determinant[1:(rep * 3),], determinant[(6 * rep + 1):(7 * rep),])
        coverage_data <- data.frame(
                            range = checkpoints,
                            coverage_ABM = sapply(1:length(checkpoints), function(i) mean(coverage[((i - 1) * rep + 1):(i * rep), 1])),
                            coverage_SM = sapply(1:length(checkpoints), function(i) mean(coverage[((i - 1) * rep + 1):(i * rep), 2])),
                            coverage_RBM = sapply(1:length(checkpoints), function(i) mean(coverage[((i - 1) * rep + 1):(i * rep), 3])),
                            coverage_TRUE = sapply(1:length(checkpoints), function(i) mean(coverage[((i - 1) * rep + 1):(i * rep), 4]))
                        )
        foo <- paste("MVG_coverage", "rho", rho, "m", m, ".pdf", sep = "_")
        foo <- paste(foo, sep = "/")
        pdf(foo, height = 11, width = 10)
        grid.table(coverage_data)
        dev.off()

    }


}


A <- 1
B <- 1
n <- 1e4
rep <- 1e3
m_ <- c(5, 10)
plot_step <- 50
checkpoints <- c(1e2, 5e2, 1e3, 1e4)
p <- 2


type = 'coverage'
for (i in 1:length(m_)) {
    m <- m_[i]

    rho <- 0.999
    make_plots(A, B, rho, m, n, rep, plot_step = plot_step, type)

    rho <- 0.5
    make_plots(A, B, rho, m, n, rep, plot_step = plot_step, type)

}





######################################
##### Rosenbrock
######################################

## CONTOUR PLOT OF ROSENBROCK DENSITY

rosenbrock_density <- function(x1, x2) {
    z = exp(-5 * (x2 - x1 ^ 2) ^ 2 - (1 / 20) * (1 - x1) ^ 2)
}

make_contour_plot <- function(m, edge1, mu, A, B) {

    z = matrix(, nrow = m, ncol = m)
    x1 <- seq(-edge1, edge1, length.out = m)
    x2 <- seq(0, 120, length.out = m)
    for (i in 1:m) {
        for (j in 1:m) {
            z[i, j] = rosenbrock_density(x1[i], x2[j])
            #z[i, j] = dnorm(x1[i], mean = mu, sd = sqrt(1 / (2 * A))) * dnorm(x2[j], mean = x1[i] ^ 2, sd = sqrt(1 / (2 * B)))
        }
    }
    plot_name <- "rosenbrock_contour.pdf"
    pdf(plot_name, height = 6, width = 6)
    persp3D(x = x1, y = x2, z = z, theta = 110, phi = 40, axes = TRUE, nticks = 5, scale = 2, box = TRUE, ticktype = "detailed", zlab = "density")
    dev.off()

}

m = 4e2
edge1 = 8
A = 1 / 20
B = 5
mu = 1
make_contour_plot(m, edge1, mu, A, B)

## TRACE PLOTS OF MARKOV CHAIN WITH ROSENBROCK AS TARGET

# log unnormalized density
lufac <- function(A, B, mu) function(X) {
    r <- -A * (X[1] - mu) ^ 2 - B * (X[2] - X[1] ^ 2) ^ 2
    return(r)
}

make_trace_plots <- function(m, n, edge1, edge2, p) {

    init <- cbind(seq(-edge1, edge1, length.out = m), sapply(seq(-edge2, edge2, length.out = m), function(t) rnorm(1, mean = t ^ 2, sd = sqrt(1 / (2 * B)))))
    lupost = lufac(A, B, mu)
    X <- array(, dim = c(m, n, p))
    for (i in 1:m) {
        X[i,,] <- metrop(lupost, initial = init[i,], nbatch = n, scale = c(1 / 2, 3 / 2))$batch
    }

    plot_name <- "rosenbrock_trace_plots.pdf"
    pdf(plot_name, height = 6, width = 6)
    yrange = range(X[1,, 2], X[2,, 2])
    plot.ts(X[1,, 2], col = alpha("blue", 1), type = 'l', lty = 1, ylim = c(yrange[1] / 1.05, yrange[2] * 1.05), ylab = 'Second component', xlab = 'Chain Length', cex.axis = 1.4, cex.lab = 1.4)
    lines(X[2,, 2], col = alpha("red", 0.8), lty = 2)

    op <- par(cex = 1.2)
    legend("topright", legend = c("chain 1", "chain 2"), col = c("blue", "red"), lty = c(1, 2), bty = 'n')
    dev.off()
}

m = 2
n = 1e5
edge1 = 8
edge2 <- 8
A = 1 / 20
B = 5
mu = 1
p = 2
make_trace_plots(m, n, edge1, edge2, p)

## PLOT FUNCTION FOR ESS RUNNING PLOTS
make_plots <- function(m, n, rep, plot_step = 100) {

    filename = paste("rosenbrock_running_ESS", edge1, edge2, "m", m, "n", n, ".Rdata", sep = "_")
    load(filename)
    starting_buffer <- 0
    n <- 2e4
    ind <- (1:(floor((n - starting_buffer) / plot_step)))

    x <- (0:(floor((n - starting_buffer) / plot_step) - 1)) * plot_step + starting_buffer
    plot_name <- paste("rosenbrock_running_ESS", edge1, edge2, "m", m, "n", n, ".pdf", sep = "_")
    pdf(plot_name, height = 6, width = 6)
    par(mar = c(5, 5, 4, 1) + .1)
    plot(x, colMeans(plot_mat_abm[, ind]), type = "l", col = "green", lty = 1, lwd = 2, ylim = range(colMeans(plot_mat_abm[, ind]), colMeans(plot_mat_rbm[, ind])), ylab = expression(hat(ESS) / mn), xlab = "Chain Length", cex.axis = 1.4, cex.lab = 1.4)

    se_abm = sqrt(apply(plot_mat_abm[, ind], 2, var) / dim(plot_mat_abm)[1])
    segments(x0 = x, y0 = colMeans(plot_mat_abm[, ind]) - 1.96 * se_abm, x1 = x, y1 = colMeans(plot_mat_abm[, ind]) + 1.96 * se_abm, col = "green")

    lines(x, colMeans(plot_mat_rbm[, ind]), type = "l", col = "blue", lty = 6, lwd = 2)
    se_rbm = sqrt(apply(plot_mat_rbm[, ind], 2, var) / dim(plot_mat_rbm)[1])
    segments(x0 = x, y0 = colMeans(plot_mat_rbm[, ind]) - 1.96 * se_rbm, x1 = x, y1 = colMeans(plot_mat_rbm[, ind]) + 1.96 * se_rbm, col = "blue")

    op <- par(cex = 1.2)
    legend("topright", legend = c("ABM", "RBM"), col = c("green","blue"), lty = c(1, 6), bty = 'n')
    dev.off()

}
edge1 = 8
edge2 = 8
rep <- 1e2 
m_ <- c(5, 10)
n <- 1e5

plot_step <- 500

for (iter1 in 1:length(m_)) {

    m <- m_[iter1]
    make_plots(m, n, rep, plot_step)
    
}


## MAKE TABLES
make_plots <- function(rep, m, n) {

    filename = paste("rosenbrock_CDE", edge1, edge2, "m", m, "n", n, ".Rdata", sep = "_")
    load(filename)

    # Tables:
    coverage = rbind(coverage[1:(2*rep),], coverage[(3*rep + 1):(4*rep),], coverage[(5*rep + 1):(6*rep),])
    coverage_data <- data.frame(
            range = checkpoints,
            coverage_ABM = sapply(1:length(checkpoints), function(i) mean(coverage[((i - 1) * rep + 1):(i * rep), 1])),
            coverage_RBM = sapply(1:length(checkpoints), function(i) mean(coverage[((i - 1) * rep + 1):(i * rep), 3]))
            )
    foo <- paste("rosennbrock_coverage", edge1, edge2, "n", n, "m", m, ".pdf", sep = "_")
    pdf(foo, height = 11, width = 10)
    grid.table(coverage_data)
    dev.off()


}

edge1 = 8
edge2 = 8
rep <- 1e3 
m_ <- c(5, 10)
n <- 1e5
plot_step <- 500
p <- 2
checkpoints <- c(5e3, 1e4, 5e4, n)


for (i in 1:length(m_)) {
    m <- m_[i]
    make_plots(rep, m, n)
}


######################################
##### Multinomial
######################################

## ACF PLOTS FOR MULTINOMIAL
data(Nethvote)
make_acf_plots <- function(m, n, p) {

    MLE <- MCMCmnl(vote ~ choicevar(distD66, "sqdist", "D66") +
    choicevar(distPvdA, "sqdist", "PvdA") +
    choicevar(distVVD, "sqdist", "VVD") +
    choicevar(distCDA, "sqdist", "CDA") +
    relig + class + income + educ + age + urban,
    baseline = "D66", mcmc.method = "RWM", B0 = 0,
    verbose = FALSE, mcmc = 1, thin = 1, tune = 0.5, burnin = 0, seed = 1,
    data = Nethvote)

    starting <- matrix(, nrow = m, ncol = p)
    for (i in 1:p) {
        starting[, i] <- seq(MLE[1, i] - 2, MLE[1, i] + 2, length = m)
    }


    Y <- MCMCmnl(vote ~ choicevar(distD66, "sqdist", "D66") +
    choicevar(distPvdA, "sqdist", "PvdA") +
    choicevar(distVVD, "sqdist", "VVD") +
    choicevar(distCDA, "sqdist", "CDA") +
    relig + class + income + educ + age + urban,
    baseline = "D66", beta.start = starting[1,], mcmc.method = "RWM", B0 = 0,
    verbose = 0, mcmc = n, thin = 1, tune = .5, burnin = 0, seed = 1,
    data = Nethvote)


    plot_name <- "multinomial_acf_plot_4.pdf"
    pdf(plot_name, height = 4, width = 6)
    acf(Y[, 4], main = " ", lag.max = 100, cex.axis = 1.4, cex.lab = 1.4)
    dev.off()

    plot_name <- "multinomial_acf_plot_18.pdf"
    pdf(plot_name, height = 4, width = 6)
    acf(Y[, 18], main = " ", lag.max = 100, cex.axis = 1.4, cex.lab = 1.4)
    dev.off()

}

m <- 2
n <- 2e4
p <- 22
make_acf_plots(m, n, p)

## PLOT FUNCTION FOR ESS RUNNING PLOTS
make_plots <- function(m, n, rep, plot_step = 100) {

    filename = paste("multinomial_running_sqroot_2", "m", m, "n", n, "plot_step", plot_step, ".Rdata", sep = "_")
    load(filename)
    starting_buffer <- 500

    ind <- (1:(1 + floor((n - starting_buffer) / plot_step)))

    x <- (0:(floor((n - starting_buffer) / plot_step))) * plot_step + starting_buffer

    plot_name <- paste("multinomial_running_sqroot_2", "m", m, "n", n, "plot_step", plot_step, ".pdf", sep = "_")
    pdf( plot_name, height = 6, width = 6)
    par(mar = c(5, 5, 4, 1) + .1)
    plot(x, colMeans(plot_mat_abm[, ind]), type = "l", col = "green", lty = 1, lwd = 2, ylim = range(colMeans(plot_mat_abm[, ind]), colMeans(plot_mat_rbm[, ind])), ylab = expression(hat(ESS)/mn), xlab = "Chain Length", cex.axis = 1.4, cex.lab = 1.4)
    se_abm = sqrt(apply(plot_mat_abm[, ind], 2, var) / dim(plot_mat_abm)[1])
    segments(x0 = x, y0 = colMeans(plot_mat_abm[, ind]) - 1.96 * se_abm, x1 = x, y1 = colMeans(plot_mat_abm[, ind]) + 1.96 * se_abm, col = "green")

    lines(x, colMeans(plot_mat_rbm[, ind]), type = "l", col = "blue", lty = 6, lwd = 2)
    se_rbm = sqrt(apply(plot_mat_rbm[, ind], 2, var) / dim(plot_mat_rbm)[1])
    segments(x0 = x, y0 = colMeans(plot_mat_rbm[, ind]) - 1.96 * se_rbm, x1 = x, y1 = colMeans(plot_mat_rbm[, ind]) + 1.96 * se_rbm, col = "blue")

    op <- par(cex = 1.2)
    legend("topright", legend = c("ABM", "RBM"), col = c("green", "blue"), lty = c(1, 6), bty = 'n')
    dev.off()

}

rep <- 1e3 
m <- 2
n <- 1e4
plot_step <- 100

make_plots(m, n, rep, plot_step)


