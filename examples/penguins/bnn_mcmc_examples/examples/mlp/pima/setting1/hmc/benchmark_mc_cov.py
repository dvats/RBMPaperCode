# %% Load packages

import numpy as np

from kanga.chains import ChainArrays

from bnn_mcmc_examples.examples.mlp.pima.setting1.constants import diagnostic_iter_thres, num_chains
from bnn_mcmc_examples.examples.mlp.pima.setting1.hmc.constants import sampler_output_path, sampler_output_run_paths

# %% Load chain arrays

chain_arrays = ChainArrays.from_file(sampler_output_run_paths, keys=['sample'])

# %% Drop burn-in samples

chain_arrays.vals['sample'] = chain_arrays.vals['sample'][:, diagnostic_iter_thres:, :]

# %% Compute Monte Carlo covariance matrices

mc_cov_mats = chain_arrays.mc_cov()

# %% Save Monte Carlo covariance matrices

for i in range(num_chains):
    np.savetxt(sampler_output_run_paths[i].joinpath('mc_cov.csv'), mc_cov_mats[i], delimiter=',')

# %% Save mean of Monte Carlo covariance matrices

np.savetxt(sampler_output_path.joinpath('mean_mc_cov.csv'), mc_cov_mats.mean(0), delimiter=',')

# %% Compute Monte Carlo correlation matrices

mc_cor_mats = chain_arrays.mc_cor(mc_cov_mat=mc_cov_mats)

# %% Save Monte Carlo correlation matrices

for i in range(num_chains):
    np.savetxt(sampler_output_run_paths[i].joinpath('mc_cor.csv'), mc_cor_mats[i], delimiter=',')

# %% Save mean of Monte Carlo correlation matrices

np.savetxt(sampler_output_path.joinpath('mean_mc_cor.csv'), mc_cor_mats.mean(0), delimiter=',')
