# %% Import packages

from bnn_mcmc_examples.examples.mlp.hawks.constants import num_chains, num_epochs, num_burnin_epochs, verbose, verbose_step
from bnn_mcmc_examples.examples.mlp.hawks.metropolis_hastings.constants import sampler_output_path
from bnn_mcmc_examples.examples.mlp.hawks.metropolis_hastings.sampler import sampler

# %% Benchmark Metropolis-Hastings sampler

sampler.benchmark(
    num_chains=num_chains,
    num_epochs=num_epochs,
    num_burnin_epochs=num_burnin_epochs,
    path=sampler_output_path,
    check_conditions=lambda chain, runtime : 0.10 <= chain.acceptance_rate() <= 0.70,
    verbose=verbose,
    verbose_step=verbose_step,
    print_acceptance=True,
    print_runtime=True
)
