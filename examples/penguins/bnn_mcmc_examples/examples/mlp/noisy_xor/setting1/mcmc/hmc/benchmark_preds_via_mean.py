# %% Load packages

import numpy as np
import torch

from eeyore.chains import ChainLists

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.constants import dtype
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.constants import num_chains, pred_iter_thres
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.datascanners import test_dataloader
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.hmc.constants import sampler_output_run_paths
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.model import model

# %% Load chain lists

chain_lists = ChainLists.from_file(sampler_output_run_paths, keys=['sample'], dtype=dtype)

# %% Drop burn-in samples

for i in range(num_chains):
    chain_lists.vals['sample'][i] = chain_lists.vals['sample'][i][pred_iter_thres:]

# %% Compute chain means

means = chain_lists.mean()

# %% Make and save predictions

for k in range(num_chains):
    test_pred_probs = np.empty([len(test_dataloader)])

    for i, (x, _) in enumerate(test_dataloader):
        integral, _ = model.predictive_posterior([means[k, :]], x, torch.tensor([[1.]], dtype=dtype))
        test_pred_probs[i] = integral.item()

    test_preds = test_pred_probs > 0.5

    np.savetxt(sampler_output_run_paths[k].joinpath('preds_via_mean.txt'), test_preds, fmt='%d')
