# %% Import packages

import kanga.plots as ps

from kanga.chains import ChainArray

from bnn_mcmc_examples.examples.mlp.mnist.constants import diagnostic_iter_thres
from bnn_mcmc_examples.examples.mlp.mnist.hmc.constants import sampler_output_pilot_path

# %% Load chain array

chain_array = ChainArray.from_file(keys=['sample', 'accepted'], path=sampler_output_pilot_path)

# %% Drop burn-in samples

chain_array.vals['sample'] = chain_array.vals['sample'][diagnostic_iter_thres:, :]
chain_array.vals['accepted'] = chain_array.vals['accepted'][diagnostic_iter_thres:]

# %% Set subset of parameters for visual summaries

par_subset = [1, 3000, 8000]
par_subset[:] = [i-1 for i in par_subset]

# %% Plot traces of simulated chain

for i in par_subset:
    ps.trace(
        chain_array.get_param(i),
        title=r'Traceplot of $\theta_{{{}}}$'.format(i+1),
        xlabel='Iteration',
        ylabel='Parameter value'
    )

# %% Plot running means of simulated chain

for i in par_subset:
    ps.running_mean(
        chain_array.get_param(i),
        title=r'Running mean plot of parameter $\theta_{{{}}}$'.format(i+1),
        xlabel='Iteration',
        ylabel='Running mean'
    )

# %% Plot histograms of marginals of simulated chain

for i in par_subset:
    ps.hist(
        chain_array.get_param(i),
        bins=30,
        density=True,
        title=r'Histogram of parameter $\theta_{{{}}}$'.format(i+1),
        xlabel='Parameter value',
        ylabel='Parameter relative frequency'
    )
