# %% Load packages

import numpy as np

from eeyore.chains import ChainLists

from bnn_mcmc_examples.examples.mlp.pima.setting2.constants import dtype, num_chains, pred_iter_thres
from bnn_mcmc_examples.examples.mlp.pima.setting2.dataloaders import test_dataloader
from bnn_mcmc_examples.examples.mlp.pima.setting2.metropolis_hastings.constants import sampler_output_run_paths
from bnn_mcmc_examples.examples.mlp.pima.setting2.model import model

# %% Load chain lists

chain_lists = ChainLists.from_file(sampler_output_run_paths, keys=['sample'], dtype=dtype)

# %% Drop burn-in samples

for i in range(num_chains):
    chain_lists.vals['sample'][i] = chain_lists.vals['sample'][i][pred_iter_thres:]

# %% Load test data and labels

test_data, test_labels = next(iter(test_dataloader))

# %% Compute chain means

means = chain_lists.mean()

# %% Make and save predictions

for i in range(num_chains):
    # Initialize model parameters
    model.set_params(means[i, :].clone().detach())

    # Compute test logits
    test_logits = model(test_data)

    # Make test predictions
    test_preds = test_logits.squeeze() > 0.5

    # Save predictions
    np.savetxt(sampler_output_run_paths[i].joinpath('preds_via_mean.txt'), test_preds, fmt='%d', delimiter=',')
