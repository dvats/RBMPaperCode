# %% Import packages

from bnn_mcmc_examples.datasets import load_xydataset_from_file
from bnn_mcmc_examples.datasets.noisy_xor.data1.constants import test_data_path
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.constants import dtype

# %% Load test dataloader with batch size of 1

test_dataset, test_dataloader = load_xydataset_from_file(test_data_path, dtype=dtype, batch_size=1)
