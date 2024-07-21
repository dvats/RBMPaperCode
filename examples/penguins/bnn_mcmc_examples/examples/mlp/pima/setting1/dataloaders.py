# %% Import packages

from torch.utils.data import DataLoader

from bnn_mcmc_examples.examples.mlp.pima.setting1.datasets import training_dataset, test_dataset

# %% Create training dataloader

training_dataloader = DataLoader(training_dataset, batch_size=len(training_dataset), shuffle=False)

# %% Create test dataloader

test_dataloader = DataLoader(test_dataset, batch_size=len(test_dataset), shuffle=False)
