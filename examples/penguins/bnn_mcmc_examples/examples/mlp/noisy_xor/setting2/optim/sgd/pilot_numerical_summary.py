# %% Import packages

import numpy as np
import torch

from sklearn.metrics import accuracy_score

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.constants import dtype
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.model import model
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.optim.dataloaders import test_dataloader
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.optim.sgd.constants import optimizer_output_pilot_path

# %% Load solution

with open(optimizer_output_pilot_path.joinpath('solution.csv'), 'r') as file:
    solution = np.loadtxt(file)

# %% Set model parameters

model.set_params(torch.tensor(solution, dtype=dtype))

# %% Load test data and labels

test_data, test_labels = next(iter(test_dataloader))

# %% Compute test logits

test_logits = model(test_data)

# Make test predictions

test_preds = test_logits.squeeze() > 0.5

# %% Compute test accuracy

print(accuracy_score(test_preds, test_labels.squeeze()))
# (preds == labels.squeeze()).sum()/len(labels.squeeze())

# %% Show parameter values

print(model.get_params())
