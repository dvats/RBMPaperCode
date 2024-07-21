# %% Import packages

import copy
import numpy as np
import torch

from eeyore.kernels import NormalKernel
from eeyore.samplers import PowerPosteriorSampler

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.dataloaders import training_dataloader
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.model import model

# %% Setup proposal variance and proposal kernel for Metropolis-Hastings per-chain samplers

proposal_scale = np.sqrt(0.02)

kernel = NormalKernel(
    torch.zeros(model.num_params(), dtype=model.dtype),
    torch.full([model.num_params()], proposal_scale, dtype=model.dtype)
)

# %% Setup power posterior sampler

num_power_posteriors = 10

per_chain_samplers = [['MetropolisHastings', {'kernel': copy.deepcopy(kernel)}] for _ in range(num_power_posteriors)]

sampler = PowerPosteriorSampler(
    model,
    training_dataloader,
    per_chain_samplers,
    theta0=model.prior.sample(),
    temperature=[1. for _ in range(num_power_posteriors)],
    between_step=1,
    keys=['sample', 'target_val'],
    check_input=True
)
