
set.seed(10)
library(reticulate)
#install_miniconda(path = miniconda_path(), update = TRUE, force = FALSE)
use_condaenv('r-reticulate', required = TRUE)

py_install(c('torch', 'sklearn'), pip = TRUE)
conda_install('r-reticulate', packages = c("eeyore", "kanga"))

source_python('reticulate.py')

########################################################################
## Aim : Calculate Variance using average batch means, 
##       replicated batch means and stupid method for 2d Rosenbrock density
########################################################################
set.seed(11)
library(mcmcse)


source("./../get_batch.R")
source("./../ABM.R")
source("./../naive.R")
source("./../RBM.R")
source("./../choosingBatch.R")
source("./../ESS.R")
# source('../batchSize_final.R')

# Log unnormalized posterior



running_penguins <- function(rep, m, n, plot_step) {
  
  starting_buffer = 0
  
  plot_mat_abm <- matrix(, nrow = rep, ncol = 1 + floor((n - starting_buffer) / plot_step))
  plot_mat_nm <- matrix(, nrow = rep, ncol = floor((n - starting_buffer) / plot_step))
  plot_mat_rbm <- matrix(, nrow = rep, ncol = 1 + floor((n - starting_buffer) / plot_step))
  
  # stores Markov chains
  X <- array(, dim = c(m, n, p))
  #init <- cbind(seq(-edge1, edge1, length.out = m), sapply(seq(-edge2, edge2, length.out = m), function(t) rnorm(1, mean = t ^ 2, sd = sqrt(1 / (2 * B)))))
  
  print(n + m)
  for (j in 1:rep) {
    
    if ((j %% 10) == 0)
      print(j)

    # sample Markov chains
    for (i in 1:m) {
      # MC is jumping rarely
      X[i,,] <- get_samples(n, verbose = FALSE)
      
    }
    
    # calculate estimates across plot steps
    for (i in (1:floor((n - starting_buffer) / plot_step))) {
      if ((i %% 10) == 0)
        print(paste("Rep = ", j, ", % = ", i / (floor(n / plot_step))))
      b <- chooseBatch(X[, 1:(starting_buffer + plot_step * i),], r,'smart')
      abm <- ABM(X[, 1:(starting_buffer + plot_step * i),], m, starting_buffer + plot_step * i, r, b)
      nm <- SM(X[, (starting_buffer + 1):(starting_buffer + plot_step * i),], m, plot_step * i)
      rbm <- RBM(X[, 1:(starting_buffer + plot_step * i),], m, starting_buffer + plot_step * i, r, c, b)
      
      plot_mat_abm[j, i] = (ESS(X[, 1:(starting_buffer + plot_step * i),], abm))[2]
      plot_mat_nm[j, i] = (ESS(X[, 1:(starting_buffer + plot_step * i),], nm))[2]
      plot_mat_rbm[j, i] = (ESS(X[, 1:(starting_buffer + plot_step * i),], rbm))[2]
    }
    
  }
  
  filename = paste("penguins_running_ESS_stationary_newadjsm","m", m, "n", n, ".Rdata", sep = "_")
  save(plot_mat_abm, plot_mat_nm, plot_mat_rbm, file = paste("output_files", filename, sep = "/"))
  
}


rep <- 2 #1e2 
m_ <- c(5, 10)
n_ <- c(1e3) #c(1e5)
p <- 29
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
    running_penguins(rep, m, n, plot_step)
    
  }
}
# time = 26440




