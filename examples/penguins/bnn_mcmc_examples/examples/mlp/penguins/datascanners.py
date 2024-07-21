# %% Import packages

from torch.utils.data import DataLoader

from bnn_mcmc_examples.examples.mlp.penguins.datasets import test_dataset

# %% Load test dataloader with batch size of 1

test_dataloader = DataLoader(test_dataset, batch_size=1, shuffle=False)
