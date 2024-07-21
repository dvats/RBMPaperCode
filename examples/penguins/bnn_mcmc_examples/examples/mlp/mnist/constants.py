# %% Import packages

import torch

from pathlib import Path

# %% Define constants

output_path = Path.home().joinpath('output', 'bnn_mcmc_examples', 'mlp', 'mnist')

num_features = 784
num_classes = 10

mlp_dims = [num_features, 10, 10, 10, num_classes]
mlp_bias = 4 * [True]
mlp_activations = 3 * [torch.sigmoid]
mlp_activations.append(None)

dtype = torch.float32

num_chains = 10

num_epochs = 1100 # 11000
num_burnin_epochs = 0
diagnostic_iter_thres = 100 # 1000

verbose = True
verbose_step = 1

pred_iter_thres = 100 # 1000
