# %% Import packages

import torchvision.datasets as datasets

from bnn_mcmc_examples.datasets import data_root

# %% Download training dataset

datasets.MNIST(root=data_root, train=True, download=True, transform=None)

# %% Download test dataset

datasets.MNIST(root=data_root, train=False, download=True, transform=None)
