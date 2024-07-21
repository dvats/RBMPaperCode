# %% Import packages

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.constants import (
    num_chains, num_epochs, num_burnin_epochs, verbose, verbose_step
)
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.metropolis_hastings.constants import sampler_output_path
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.metropolis_hastings.sampler import sampler

# %% Benchmark Metropolis-Hastings sampler

sampler.benchmark(
    num_chains=num_chains,
    num_epochs=num_epochs,
    num_burnin_epochs=num_burnin_epochs,
    path=sampler_output_path,
    check_conditions=lambda chain, runtime : 0.15 <= chain.acceptance_rate() <= 0.50,
    verbose=verbose,
    verbose_step=verbose_step,
    print_acceptance=True,
    print_runtime=True
)
