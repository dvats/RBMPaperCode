# %% Import packages

from torch.utils.data import DataLoader

from eeyore.datasets import XYDataset

from bnn_mcmc_examples.examples.mlp.exact_xor.constants import dtype

# %% Load dataloader

dataset = XYDataset.from_eeyore('xor', dtype=dtype)

dataloader = DataLoader(dataset, batch_size=1)
