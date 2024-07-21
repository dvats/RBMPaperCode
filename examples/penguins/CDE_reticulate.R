set.seed(10)
library(mcmcse)
source("./../ABM.R")
source("./../naive.R")
source("./../RBM.R")
source("./../get_batch.R")
source("./../ESS.R")
source("./../choosingBatch.R")

library(reticulate)
#install_miniconda(path = miniconda_path(), update = TRUE, force = FALSE)
use_condaenv('r-reticulate', required = TRUE)

py_install(c('torch', 'sklearn'), pip = TRUE)
conda_install('r-reticulate', packages = c("eeyore", "kanga"))

source_python('reticulate.py')


running_penguins <- function(rep, m, n) {
  print(n + m)

  # stores Markov chains
  X <- array(, dim = c(m, n, p))
  
  #determinant <- matrix(0, nrow = length(checkpoints) * rep, ncol = 3)
  ess_mat_concat <- matrix(, nrow = length(checkpoints) * rep, ncol = 3)
  ess_mat_avg <- matrix(, nrow = length(checkpoints) * rep, ncol = 3)
  
  ABM_estimate <- array(dim = c(length(checkpoints) * rep, p, p))
  SM_estimate <- array(dim = c(length(checkpoints) * rep, p, p))
  RBM_estimate <- array(dim = c(length(checkpoints) * rep, p, p))
  
  for (j in 1:rep) {
    
    if ((j %% (j/10)) == 0)
      print(j)
    
    # First draw of the Markov chain, 
    
    for (i in 1:m) {
      X[i,,] <- get_samples(n, TRUE)
    }
    
    # calculate covariance estimate at each checkpoint
    for (i in 1:length(checkpoints)) {
      b <- chooseBatch(X[, 1:checkpoints[i],], r, 'sqroot')
      ABM_estimate[(i - 1) * rep + j,,] <- ABM(X[, 1:checkpoints[i],], m, checkpoints[i], r, b)
      SM_estimate[(i - 1) * rep + j,,] <- SM(X[, 1:checkpoints[i],], m, checkpoints[i])
      RBM_estimate[(i - 1) * rep + j,,] <- RBM(X[, 1:checkpoints[i],], m, checkpoints[i], r, c, b)
      
      
      #determinant[(i - 1) * rep + j, 1] = det(ABM_estimate[(i - 1) * rep + j,,])
      #determinant[(i - 1) * rep + j, 2] = det(SM_estimate[(i - 1) * rep + j,,])
      #determinant[(i - 1) * rep + j, 3] = det(RBM_estimate[(i - 1) * rep + j,,])
      
      ess_abm <- ESS(X[, 1:checkpoints[i],], ABM_estimate[(i - 1) * rep + j,,])
      ess_mat_concat[(i - 1) * rep + j, 1] <- ess_abm[1]
      ess_mat_avg[(i - 1) * rep + j, 1] <- ess_abm[2]
      
      ess_sm <- ESS(X[, 1:checkpoints[i],], SM_estimate[(i - 1) * rep + j,,])
      ess_mat_concat[(i - 1) * rep + j, 2] <- ess_sm[1]
      ess_mat_avg[(i - 1) * rep + j, 2] <- ess_sm[2]
      
      ess_rbm <- ESS(X[, 1:checkpoints[i],], RBM_estimate[(i - 1) * rep + j,,])
      ess_mat_concat[(i - 1) * rep + j, 3] <- ess_rbm[1]
      ess_mat_avg[(i - 1) * rep + j, 3] <- ess_rbm[2]
      
    }
    
    
  }
  filename = paste("CDE_penguins", "m", m, "n", n, ".Rdata", sep = "_")
  save(ABM_estimate, SM_estimate, RBM_estimate, ess_mat_concat, ess_mat_avg, file = paste("output_files", filename, sep = "/"))
  
}

rep <- 1e2 
m_ <- c(5)
n <- 1e4
alpha <- 0.05
p <- 29
plot_step <- 500
# checkpoints <- c(5e3, 1e4, 3e4, 5e4, 7e4, n)
checkpoints <- c(5e2, 1e3, 2e3, 5e3, 7e3, n)
r <- 3 # Lugsail hyperparameters
c <- 0.5 # Lugsail hyperparameters

start.time <- Sys.time()

for (i in 1:length(m_)) {
  m <- m_[i]
  running_penguins(rep, m, n)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)






