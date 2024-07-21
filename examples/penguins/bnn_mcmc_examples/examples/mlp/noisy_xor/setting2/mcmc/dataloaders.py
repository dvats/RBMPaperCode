# %% Import packages

from bnn_mcmc_examples.datasets import load_xydataset_from_file
from bnn_mcmc_examples.datasets.noisy_xor.data2.constants import test_data_path, training_data_path
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.constants import dtype

# %% Load training dataloader

training_dataset, training_dataloader = load_xydataset_from_file(training_data_path, dtype=dtype)

# %% Load test dataloader

test_dataset, test_dataloader = load_xydataset_from_file(test_data_path, dtype=dtype)
