# %% Load packages

import numpy as np

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.constants import output_path
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.constants import pred_interval_x1, pred_interval_x2

# %% Generate ground truth

pred_interval_y = np.empty([len(pred_interval_x1), len(pred_interval_x2)])

for i in range(len(pred_interval_x1)):
    for j in range(len(pred_interval_x2)):
        if pred_interval_x1[i] < 0.5:
            pred_interval_y[i, j] = 0 if (pred_interval_x2[j] < 0.5) else 1
        else:
            pred_interval_y[i, j] = 1 if (pred_interval_x2[j] < 0.5) else 0

# %% Save ground truth

np.savetxt(output_path.joinpath('ground_truth.csv'), pred_interval_y, delimiter=',')
