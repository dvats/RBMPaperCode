# %% Import packages

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.constants import (
    num_chains, num_epochs, num_burnin_epochs, verbose, verbose_step
)
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.power_posteriors.constants import sampler_output_path
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.power_posteriors.sampler import sampler

# %% Benchmark power posterior sampler

sampler.benchmark(
    num_chains=num_chains,
    num_epochs=num_epochs,
    num_burnin_epochs=num_burnin_epochs,
    path=sampler_output_path,
    check_conditions=None,
    verbose=verbose,
    verbose_step=verbose_step,
    print_acceptance=False,
    print_runtime=True
)
