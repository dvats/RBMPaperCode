# %% Import packages

from eeyore.samplers import HMC

from bnn_mcmc_examples.examples.mlp.pima.setting1.dataloaders import training_dataloader
from bnn_mcmc_examples.examples.mlp.pima.setting1.model import model

# %% Setup HMC sampler

sampler = HMC(model, theta0=model.prior.sample(), dataloader=training_dataloader, step=0.125, num_steps=6)
