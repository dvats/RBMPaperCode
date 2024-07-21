
set.seed(10)

########################################################################
## Aim : Calculate Variance using average batch means, 
##       replicated batch means and stupid method for 2d Rosenbrock density
########################################################################
set.seed(11)
library(mcmcse)

options("scipen"=10)
source("./../get_batch.R")
source("./../ABM.R")
source("./../naive.R")
source("./../RBM.R")
source("./../choosingBatch.R")
source("./../ESS.R")
# source('../batchSize_final.R')

# Log unnormalized posterior



running_penguins <- function(rep, m, n, plot_step) {
  
  starting_buffer = 1e3
  
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
      filename = paste('Mchains/m_',m,'_n_',n,'_/chain_m__',j-1,'_',i-1,'.csv',sep = "")
      #filename = paste('Mchains/m_',m,"_n_",n,'_/chain_m_',"_",j-1,"_",i-1,'.csv',sep = "")
      X[i,,] <- as.matrix(read.csv(filename, header = FALSE)) #get_samples(n, verbose = FALSE)
      
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
    filename = paste("penguins_running_ESS_stationary_newadjsm","m", m, "n", n, ".Rdata", sep = "_")
    save(plot_mat_abm, plot_mat_nm, plot_mat_rbm, file = paste("output_files", filename, sep = "/"))
  }
  
  filename = paste("penguins_running_ESS_stationary_newadjsm","m", m, "n", n, ".Rdata", sep = "_")
  save(plot_mat_abm, plot_mat_nm, plot_mat_rbm, file = paste("output_files", filename, sep = "/"))
  
}


rep <- 10
m_ <- c(5) #, 10)
n_ <- c(1e5) #c(1e5)
p <- 29
plot_step <- 2000
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




