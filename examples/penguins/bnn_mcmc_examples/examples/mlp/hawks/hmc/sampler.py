# %% Import packages

from eeyore.samplers import HMC

from bnn_mcmc_examples.examples.mlp.hawks.dataloaders import training_dataloader
from bnn_mcmc_examples.examples.mlp.hawks.model import model

# %% Setup HMC sampler

sampler = HMC(model, theta0=model.prior.sample(), dataloader=training_dataloader, step=0.08, num_steps=6)
