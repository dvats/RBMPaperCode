# %% Load packages

import numpy as np

from kanga.chains import ChainArrays

from bnn_mcmc_examples.examples.mlp.penguins.constants import diagnostic_iter_thres, num_chains
from bnn_mcmc_examples.examples.mlp.penguins.metropolis_hastings.constants import (
    sampler_output_path, sampler_output_run_paths
)

# %% Load chain arrays, covariance matrices and runtimes

chain_arrays = ChainArrays.from_file(sampler_output_run_paths, keys=['sample'])

mc_cov_mats = []
for i in range(num_chains):
    mc_cov_mats.append(np.loadtxt(sampler_output_run_paths[i].joinpath('mc_cov.csv'), delimiter=',', skiprows=0))
mc_cov_mats = np.stack(mc_cov_mats)

runtimes = []
for i in range(num_chains):
    with open(sampler_output_run_paths[i].joinpath('runtime.txt'), 'r') as file:
        runtimes.append(float(file.readline().rstrip()))
runtimes = np.array(runtimes)

# %% Drop burn-in samples

chain_arrays.vals['sample'] = chain_arrays.vals['sample'][:, diagnostic_iter_thres:, :]

# %% Compute multivariate rhat

rhat_val, _, _, _, _ = chain_arrays.multi_rhat(mc_cov_mat=mc_cov_mats)

# %% Save rhat_val

with open(sampler_output_path.joinpath('multi_rhat.txt'), 'w') as file:
    file.write('{}\n'.format(rhat_val))

# %% Compute multivariate ESS

ess_vals = np.array(chain_arrays.multi_ess(mc_cov_mat=mc_cov_mats))

# %% Save multivariate ESSs

for i in range(num_chains):
    with open(sampler_output_run_paths[i].joinpath('multi_ess.txt'), 'w') as file:
        file.write('{}\n'.format(ess_vals[i]))

# %% Save mean of multivariate ESSs

with open(sampler_output_path.joinpath('mean_multi_ess.txt'), 'w') as file:
    file.write('{}\n'.format(ess_vals.mean()))

# %% Save mean runtime

with open(sampler_output_path.joinpath('mean_runtime.txt'), 'w') as file:
    file.write('{}\n'.format(runtimes.mean()))
