# %% Load packages

import numpy as np
import torch

from sklearn.metrics import accuracy_score

from bnn_mcmc_examples.examples.mlp.mnist.constants import num_chains
from bnn_mcmc_examples.examples.mlp.mnist.dataloaders import test_dataloader
from bnn_mcmc_examples.examples.mlp.mnist.prior.constants import sampler_output_path, sampler_output_run_paths

# %% Load test data and labels

_, test_labels = next(iter(test_dataloader))

# %% Compute predictive accuracies

accuracies = np.empty(num_chains)

for i in range(num_chains):
    test_preds = np.loadtxt(sampler_output_run_paths[i].joinpath('preds_via_mean.txt'), skiprows=0)

    accuracies[i] = accuracy_score(test_preds, torch.argmax(test_labels, 1))

# %% Save predictive accuracies

np.savetxt(sampler_output_path.joinpath('accuracies_via_mean.txt'), accuracies)
