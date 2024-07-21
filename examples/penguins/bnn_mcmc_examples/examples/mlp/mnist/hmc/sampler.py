# %% Import packages

from eeyore.samplers import HMC

from bnn_mcmc_examples.examples.mlp.mnist.dataloaders import training_dataloader
from bnn_mcmc_examples.examples.mlp.mnist.model import model

# %% Setup HMC sampler

sampler = HMC(model, theta0=model.prior.sample(), dataloader=training_dataloader, step=0.001, num_steps=6)
