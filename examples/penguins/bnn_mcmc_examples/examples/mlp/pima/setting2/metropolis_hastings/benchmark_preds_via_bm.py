# %% Load packages

import numpy as np

from bnn_mcmc_examples.examples.mlp.pima.setting2.constants import num_chains
from bnn_mcmc_examples.examples.mlp.pima.setting2.metropolis_hastings.constants import sampler_output_run_paths

# %% Make and save predictions

for i in range(num_chains):
    # Load test logits
    test_logits = np.loadtxt(sampler_output_run_paths[i].joinpath('pred_posterior_on_test.txt'), delimiter=',', skiprows=0)

    # Make test predictions
    test_preds = test_logits > 0.5

    # Save predictions
    np.savetxt(sampler_output_run_paths[i].joinpath('preds_via_bm.txt'), test_preds, fmt='%d', delimiter=',')
