# %% Import packages

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.constants import output_path

# %% Define optimizer-specific output directories

optimizer_output_path = output_path.joinpath('sgd')
optimizer_output_pilot_path = optimizer_output_path.joinpath('pilot_run')
optimizer_output_benchmark_path = optimizer_output_path.joinpath('benchmark_run')
