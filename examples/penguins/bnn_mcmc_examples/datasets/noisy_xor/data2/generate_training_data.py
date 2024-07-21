# %% Load packages

import numpy as np

from bnn_mcmc_examples.datasets.noisy_xor import simulate_xor
from bnn_mcmc_examples.datasets.noisy_xor.data2.constants import num_training_samples, training_data_path

# %% Create data directory if it does not exist

training_data_path.mkdir(parents=True, exist_ok=True)

# %% Simulate noisy XOR points

x, y = simulate_xor(n=num_training_samples)

# %% Save simulated noisy XOR points

np.savetxt(training_data_path.joinpath('x.csv'), x, delimiter=',', header='x1,x2', comments='')
np.savetxt(training_data_path.joinpath('y.csv'), y, fmt='%d', header='y', comments='')
