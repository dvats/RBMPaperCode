from eeyore.samplers import HMC
from eeyore.chains import ChainList

from bnn_mcmc_examples.examples.mlp.penguins.dataloaders import training_dataloader
from bnn_mcmc_examples.examples.mlp.penguins.model import model
from bnn_mcmc_examples.examples.mlp.penguins.constants import num_burnin_epochs, num_epochs, verbose, verbose_step

def get_samples(no_samples, verbose = False):
	no_samples = int(no_samples)
	sampler = HMC(model, theta0=model.prior.sample(), dataloader=training_dataloader, step=0.14, num_steps=6, chain = ChainList())

	sampler.run(num_epochs=no_samples, num_burnin_epochs=num_burnin_epochs, verbose=verbose, verbose_step=verbose_step)

	chain = sampler.get_chain()
	samples = chain.get_samples()
	samples = samples.numpy()
	return(samples)