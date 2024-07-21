# %% Load packages

import numpy as np
import torch

from eeyore.chains import ChainLists

from bnn_mcmc_examples.examples.mlp.hawks.constants import dtype, num_chains, num_classes, pred_iter_thres
from bnn_mcmc_examples.examples.mlp.hawks.datascanners import test_dataloader
from bnn_mcmc_examples.examples.mlp.hawks.metropolis_hastings.constants import sampler_output_run_paths
from bnn_mcmc_examples.examples.mlp.hawks.model import model

# %% Load chain lists

chain_lists = ChainLists.from_file(sampler_output_run_paths, keys=['sample'], dtype=dtype)

# %% Drop burn-in samples

for i in range(num_chains):
    chain_lists.vals['sample'][i] = chain_lists.vals['sample'][i][pred_iter_thres:]

# %% Compute chain means

means = chain_lists.mean()

# %% Make and save predictions

for k in range(num_chains):
    test_pred_probs = np.empty([len(test_dataloader), num_classes])

    for i, (x, _) in enumerate(test_dataloader):
        for j in range(num_classes):
            y = torch.zeros([1, num_classes], dtype=dtype)
            y[0, j] = 1.
            integral, _ = model.predictive_posterior([means[k, :]], x, y)
            test_pred_probs[i, j] = integral.item()

    test_preds = np.argmax(test_pred_probs, axis=1)

    np.savetxt(sampler_output_run_paths[k].joinpath('preds_via_mean.txt'), test_preds, fmt='%d')
