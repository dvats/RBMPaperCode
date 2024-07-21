# %% Load packages

import numpy as np

from bnn_mcmc_examples.examples.mlp.exact_xor.constants import num_chains
from bnn_mcmc_examples.examples.mlp.exact_xor.metropolis_hastings.constants import sampler_output_run_paths

# %% Make and save predictions

for i in range(num_chains):
    test_pred_probs = np.loadtxt(
        sampler_output_run_paths[i].joinpath('pred_posterior_on_test.csv'), delimiter=',', skiprows=0
    )

    test_preds = np.argmax(test_pred_probs, axis=1)

    np.savetxt(sampler_output_run_paths[i].joinpath('preds_via_bm.txt'), test_preds, fmt='%d')
