# %% Import packages

import json
import matplotlib.pyplot as plt

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.optim.sgd.constants import optimizer_output_pilot_path

# %% Load summary

with open(optimizer_output_pilot_path.joinpath('summary.json'), 'r') as file:
    summary = json.load(file)

# %% Plot training loss

plt.plot(summary['loss_vals'])

# %% Plot training accuracy

plt.plot(summary['metric_vals'])

# %% Plot running mean of training accuracy

plt.plot(summary['metric_mean_vals'])
