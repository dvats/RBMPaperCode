# %% Import packages

import numpy as np
import torch

from eeyore.kernels import NormalKernel
from eeyore.samplers import MetropolisHastings

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.mcmc.dataloaders import training_dataloader
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.model import model

# %% Setup proposal variance and proposal kernel for Metropolis-Hastings sampler

proposal_scale = np.sqrt(0.02)

kernel = NormalKernel(
    torch.zeros(model.num_params(), dtype=model.dtype),
    torch.full([model.num_params()], proposal_scale, dtype=model.dtype)
)

# %% Setup Metropolis-Hastings sampler

sampler = MetropolisHastings(model, theta0=model.prior.sample(), dataloader=training_dataloader, kernel=kernel)
