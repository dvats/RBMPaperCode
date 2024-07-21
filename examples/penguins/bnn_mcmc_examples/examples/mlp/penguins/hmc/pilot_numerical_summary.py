# %% Import packages

from kanga.chains import ChainArray

from bnn_mcmc_examples.examples.mlp.penguins.constants import diagnostic_iter_thres
from bnn_mcmc_examples.examples.mlp.penguins.hmc.constants import sampler_output_pilot_path

# %% Load chain array

chain_array = ChainArray.from_file(keys=['sample', 'accepted'], path=sampler_output_pilot_path)

# %% Drop burn-in samples

chain_array.vals['sample'] = chain_array.vals['sample'][diagnostic_iter_thres:, :]
chain_array.vals['accepted'] = chain_array.vals['accepted'][diagnostic_iter_thres:]

# %% Compute acceptance rate

print('Acceptance rate: {}'.format(chain_array.acceptance_rate()))

# %% Compute Monte Carlo mean

print('Monte Carlo mean: {}'.format(chain_array.mean()))

# %% Compute Monte Carlo covariance

mc_cov_mat = chain_array.mc_cov()

# %% Compute Monte Carlo standard error

print('Monte Carlo standard error: {}'.format(chain_array.mc_se(mc_cov_mat=mc_cov_mat)))

# %% Compute multivariate ESS

print('Multivariate ESS: {}'.format(chain_array.multi_ess(mc_cov_mat=mc_cov_mat)))
