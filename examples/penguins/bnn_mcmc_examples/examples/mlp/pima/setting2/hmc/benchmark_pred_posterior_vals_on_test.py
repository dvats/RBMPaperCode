# %% Load packages

import numpy as np
import torch

from eeyore.chains import ChainLists

from bnn_mcmc_examples.examples.mlp.pima.setting2.constants import dtype, num_chains, pred_iter_thres
from bnn_mcmc_examples.examples.mlp.pima.setting2.datascanners import test_dataloader
from bnn_mcmc_examples.examples.mlp.pima.setting2.hmc.constants import sampler_output_run_paths
from bnn_mcmc_examples.examples.mlp.pima.setting2.model import model

# %% Load chain lists

chain_lists = ChainLists.from_file(sampler_output_run_paths, keys=['sample'], dtype=dtype)

# %% Drop burn-in samples

for i in range(num_chains):
    chain_lists.vals['sample'][i] = chain_lists.vals['sample'][i][pred_iter_thres:]

# %% Compute and save predictive posteriors

pred_posterior = np.empty([num_chains, len(test_dataloader)])

verbose_msg = 'Evaluating predictive posterior based on chain {:' \
    + str(len(str(num_chains))) \
    + '} out of ' \
    + str(num_chains) \
    + ' at test point {:' \
    + str(len(str(len(test_dataloader)))) \
    + '} out of ' \
    + str(len(test_dataloader)) \
    + '...'

for k in range(num_chains):
    for i, (x, _) in enumerate(test_dataloader):
        print(verbose_msg.format(k+1, i+1))

        pred_posterior[k, i] = model.predictive_posterior(
            chain_lists.vals['sample'][k], x.squeeze(), torch.tensor([1.], dtype=dtype)
        )

    np.savetxt(sampler_output_run_paths[k].joinpath('pred_posterior_on_test.txt'), pred_posterior[k], delimiter=',')
