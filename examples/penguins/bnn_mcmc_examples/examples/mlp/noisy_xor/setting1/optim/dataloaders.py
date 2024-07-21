# %% Import packages

from bnn_mcmc_examples.datasets import load_xydataset_from_file
from bnn_mcmc_examples.datasets.noisy_xor.data1.constants import test_data_path, training_data_path
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.constants import dtype
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.optim.constants import batch_size, shuffle

# %% Load training dataloader

training_dataset, training_dataloader = load_xydataset_from_file(
    training_data_path, dtype=dtype, batch_size=batch_size, shuffle=shuffle
)

# %% Load test dataloader

test_dataset, test_dataloader = load_xydataset_from_file(test_data_path, dtype=dtype)
