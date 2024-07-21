# %% Load packages

import numpy as np
import torch

from eeyore.chains import ChainLists

from bnn_mcmc_examples.examples.mlp.hawks.constants import dtype, num_chains, num_classes, pred_iter_thres
from bnn_mcmc_examples.examples.mlp.hawks.datascanners import test_dataloader
from bnn_mcmc_examples.examples.mlp.hawks.model import model
from bnn_mcmc_examples.examples.mlp.hawks.prior.constants import sampler_output_run_paths

# %% Load chain lists

chain_lists = ChainLists.from_file(sampler_output_run_paths, keys=['sample'], dtype=dtype)

# %% Drop burn-in samples

for i in range(num_chains):
    chain_lists.vals['sample'][i] = chain_lists.vals['sample'][i][pred_iter_thres:]

# %% Compute and save predictive priors

verbose_msg = 'Evaluating predictive prior based on chain {:' \
    + str(len(str(num_chains))) \
    + '} out of ' \
    + str(num_chains) \
    + ' at test point {:' \
    + str(len(str(len(test_dataloader)))) \
    + '} out of ' \
    + str(len(test_dataloader)) \
    + '...'

for k in range(num_chains):
    test_pred_probs = np.empty([len(test_dataloader), num_classes])
    nums_dropped_samples = np.empty([len(test_dataloader), num_classes], dtype=np.int64)

    for i, (x, _) in enumerate(test_dataloader):
        print(verbose_msg.format(k+1, i+1))

        for j in range(num_classes):
            y = torch.zeros([1, num_classes], dtype=dtype)
            y[0, j] = 1.
            integral, num_dropped_samples = model.predictive_posterior(chain_lists.vals['sample'][k], x, y)
            test_pred_probs[i, j] = integral.item()
            nums_dropped_samples[i, j] = num_dropped_samples

    np.savetxt(sampler_output_run_paths[k].joinpath('pred_prior_on_test.csv'), test_pred_probs, delimiter=',')
    np.savetxt(
        sampler_output_run_paths[k].joinpath('pred_prior_on_test_num_dropped_samples.csv'),
        nums_dropped_samples,
        fmt='%d',
        delimiter=','
    )
