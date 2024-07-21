# %% Load packages

import numpy as np
import torch

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.constants import num_chains, num_epochs
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.prior.constants import sampler_output_run_paths
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.model import model

# %% Create output directories if they do not exist

for i in range(num_chains):
    sampler_output_run_paths[i].mkdir(parents=True, exist_ok=True)

# %% Draw samples from the prior and save them

for i in range(num_chains):
    prior_samples = torch.stack([model.prior.sample() for _ in range(num_epochs)])

    np.savetxt(sampler_output_run_paths[i].joinpath('sample.csv'), prior_samples.detach().cpu().numpy(), delimiter=',')
