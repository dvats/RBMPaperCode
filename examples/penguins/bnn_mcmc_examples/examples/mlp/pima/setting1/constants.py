# %% Import packages

import torch

from pathlib import Path

# %% Define constants

output_path = Path.home().joinpath('output', 'bnn_mcmc_examples', 'mlp', 'pima', 'setting1')

num_features = 8

mlp_dims = [num_features, 2, 2, 1]
mlp_bias = 3 * [True]
mlp_activations = 3 * [torch.sigmoid]

dtype = torch.float32

num_chains = 10

num_epochs = 110000
num_burnin_epochs = 0
diagnostic_iter_thres = 10000

verbose = True
verbose_step = 1000

pred_iter_thres = 100000
