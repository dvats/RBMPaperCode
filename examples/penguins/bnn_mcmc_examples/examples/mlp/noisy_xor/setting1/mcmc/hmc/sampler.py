# %% Import packages

from eeyore.samplers import HMC

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.dataloaders import training_dataloader
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.model import model

# %% Setup HMC sampler

sampler = HMC(model, theta0=model.prior.sample(), dataloader=training_dataloader, step=0.11, num_steps=10)
