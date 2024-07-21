# %% Import packages

from bnn_mcmc_examples.examples.mlp.mnist.constants import num_chains, num_epochs, num_burnin_epochs, verbose, verbose_step
from bnn_mcmc_examples.examples.mlp.mnist.hmc.constants import sampler_output_path
from bnn_mcmc_examples.examples.mlp.mnist.hmc.sampler import sampler

# %% Benchmark HMC sampler

sampler.benchmark(
    num_chains=num_chains,
    num_epochs=num_epochs,
    num_burnin_epochs=num_burnin_epochs,
    path=sampler_output_path,
    check_conditions=lambda chain, runtime : 0.15 <= chain.acceptance_rate(),
    verbose=verbose,
    verbose_step=verbose_step,
    print_acceptance=True,
    print_runtime=True
)
