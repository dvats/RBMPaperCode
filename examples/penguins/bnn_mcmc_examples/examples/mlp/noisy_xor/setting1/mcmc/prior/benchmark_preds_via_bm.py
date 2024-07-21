# %% Load packages

import numpy as np

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.constants import num_chains
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.prior.constants import sampler_output_run_paths

# %% Make and save predictions

for i in range(num_chains):
    test_pred_probs = np.loadtxt(
        sampler_output_run_paths[i].joinpath('pred_prior_on_test.csv'), delimiter=',', skiprows=0
    )

    test_preds = np.argmax(test_pred_probs, axis=1)

    np.savetxt(sampler_output_run_paths[i].joinpath('preds_via_bm.txt'), test_preds, fmt='%d')
