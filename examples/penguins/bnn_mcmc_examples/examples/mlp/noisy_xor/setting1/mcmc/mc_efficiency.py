# %% Load packages

import pandas as pd

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.constants import output_path
from bnn_mcmc_examples.stats import mc_efficiency

# %% Define sampler-specific output directories

sampler_output_paths = output_path.joinpath('metropolis_hastings')

sampler_output_paths = [
    output_path.joinpath(name) for name in ['metropolis_hastings', 'mala', 'smmala']
]

# %% Compute Monte Carlo efficiency

efficiency = mc_efficiency(sampler_output_paths, keys=['rhat', 'ess'])

# %% Save Monte Carlo efficiency

df = pd.DataFrame(data=efficiency)

df.to_csv(output_path.joinpath('mc_efficiency.csv'), index=False)

with open(output_path.joinpath('mc_efficiency.tex'), 'w') as file:
    file.write(df.round({'rhat': 4}).to_latex(index=False))
