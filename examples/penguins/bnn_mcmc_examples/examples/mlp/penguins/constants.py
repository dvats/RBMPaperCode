# %% Import packages

import torch

from pathlib import Path

# %% Define constants

output_path = Path.home().joinpath('output', 'bnn_mcmc_examples', 'mlp', 'penguins')

num_features = 6
num_classes = 3

mlp_dims = [num_features, 2, 2, num_classes]
mlp_bias = 3 * [True]
mlp_activations = 2 * [torch.sigmoid]
mlp_activations.append(None)

dtype = torch.float32

num_chains = 10

num_epochs = 110000
num_burnin_epochs = 0
diagnostic_iter_thres = 10000

verbose = True
verbose_step = 1000

pred_iter_thres = 100000
