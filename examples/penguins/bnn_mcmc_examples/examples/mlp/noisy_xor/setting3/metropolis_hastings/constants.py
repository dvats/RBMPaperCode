# %% Import packages

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting3.constants import num_chains, output_path

# %% Define sampler-specific output directories

sampler_output_path = output_path.joinpath('metropolis_hastings')
sampler_output_pilot_path = sampler_output_path.joinpath('pilot_run')
sampler_output_run_paths = [
    sampler_output_path.joinpath('run'+str(i+1).zfill(len(str(num_chains)))) for i in range(num_chains)
]
