# %% Load packages

import numpy as np

import kanga.plots as ps

from bnn_mcmc_examples.examples.mlp.pima.setting2.constants import num_chains
from bnn_mcmc_examples.examples.mlp.pima.setting2.metropolis_hastings.constants import (
    sampler_output_path, sampler_output_run_paths
)

# %% Load correlation matrices

mc_cor_mats = []
for i in range(num_chains):
    mc_cor_mats.append(np.loadtxt(sampler_output_run_paths[i].joinpath('mc_cor.csv'), delimiter=',', skiprows=0))
mc_cor_mats = np.stack(mc_cor_mats)

mean_mc_cor_mat = np.loadtxt(sampler_output_path.joinpath('mean_mc_cor.csv'), delimiter=',', skiprows=0)

# %% Plot heat maps of correlation matrices

for i in range(num_chains):
    ps.cor_heatmap(mc_cor_mats[i], fname=sampler_output_run_paths[i].joinpath('mc_cor.png'))

ps.cor_heatmap(mean_mc_cor_mat, fname=sampler_output_path.joinpath('mean_mc_cor.png'))
