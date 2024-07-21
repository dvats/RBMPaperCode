# %% Load packages

import numpy as np
import torch

from eeyore.chains import ChainLists

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.constants import dtype
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.mcmc.constants import (
    num_chains, pred_interval_x1, pred_interval_x2, pred_iter_thres
)
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.mcmc.metropolis_hastings.constants import sampler_output_run_paths
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.model import model

# %% Load chain lists

chain_lists = ChainLists.from_file(sampler_output_run_paths, keys=['sample'], dtype=dtype)

# %% Drop burn-in samples

for i in range(num_chains):
    chain_lists.vals['sample'][i] = chain_lists.vals['sample'][i][pred_iter_thres:]

# %% Compute and save predictive posteriors

verbose_msg = 'Evaluating predictive posterior based on chain {:' \
    + str(len(str(num_chains))) \
    + '} out of ' \
    + str(num_chains) \
    + ' at grid point ({:' \
    + str(len(str(len(pred_interval_x1)))) \
    + '}, {:' \
    + str(len(str(len(pred_interval_x2)))) \
    + '}) out of (' \
    + str(len(pred_interval_x1)) \
    + ', ' \
    + str(len(pred_interval_x2)) \
    + ')...'

for k in range(num_chains):
    pred_posterior = np.empty([len(pred_interval_x1), len(pred_interval_x2)])
    nums_dropped_samples = np.empty([len(pred_interval_x1), len(pred_interval_x2)], dtype=np.int64)

    for i in range(len(pred_interval_x1)):
        for j in range(len(pred_interval_x2)):
            print(verbose_msg.format(k+1, i+1, j+1))

            integral, num_dropped_samples = model.predictive_posterior(
                chain_lists.vals['sample'][k],
                torch.tensor([[pred_interval_x1[i], pred_interval_x2[j]]], dtype=dtype),
                torch.tensor([[1.]], dtype=dtype)
            )
            pred_posterior[i, j] = integral.item()
            nums_dropped_samples[i, j] = num_dropped_samples

    np.savetxt(sampler_output_run_paths[k].joinpath('pred_posterior_on_grid.csv'), pred_posterior, delimiter=',')
    np.savetxt(
        sampler_output_run_paths[k].joinpath('pred_posterior_on_grid_num_dropped_samples.csv'),
        nums_dropped_samples,
        fmt='%d',
        delimiter=','
    )
