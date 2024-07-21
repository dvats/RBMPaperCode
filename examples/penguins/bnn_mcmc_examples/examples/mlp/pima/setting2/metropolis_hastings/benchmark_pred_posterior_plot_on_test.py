# %% Load packages

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import torch

from bnn_mcmc_examples.examples.mlp.pima.setting2.constants import num_chains
from bnn_mcmc_examples.examples.mlp.pima.setting2.dataloaders import test_dataloader
from bnn_mcmc_examples.examples.mlp.pima.setting2.metropolis_hastings.constants import sampler_output_run_paths

# %% Load test data and labels

test_data, test_labels = next(iter(test_dataloader))

# %% Plot predictive posteriors

pred_colors = {'correct': '#bcbd22', 'wrong': '#d62728'}
# '#bcbd22': rio grande, similar to yellow green
# ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

patch_list = []
for key in pred_colors:
    patch_list.append(mpatches.Patch(color=pred_colors[key], label=key))

legend_patches = [mpatches.Patch(color=pred_colors[key], label=key.capitalize()) for key in pred_colors]

for i in range(num_chains):
    test_logits = np.loadtxt(sampler_output_run_paths[i].joinpath('pred_posterior_on_test.txt'), delimiter=',', skiprows=0)

    test_logit_dict = {
        0: np.sort(test_logits[(1 - test_labels.squeeze(-1)).to(dtype=torch.bool)]),
        1: np.sort(test_logits[test_labels.squeeze(-1).to(dtype=torch.bool)])
    }

    num_correct_0s = sum(test_logit_dict[0] < 0.5)
    num_correct_1s = sum(test_logit_dict[1] >= 0.5)
    bar_colors = num_correct_0s * [pred_colors['correct']]
    bar_colors.extend((len(test_logit_dict[0]) - num_correct_0s) * [pred_colors['wrong']])
    bar_colors.extend((len(test_logit_dict[1]) - num_correct_1s) * [pred_colors['wrong']])
    bar_colors.extend(num_correct_1s * [pred_colors['correct']])

    plt.figure(figsize=[8, 4])

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 12

    plt.vlines(
        x=range(len(test_labels)),
        ymin=0,
        ymax=np.hstack([test_logit_dict[0], test_logit_dict[1]]),
        color=bar_colors,
        linewidth=2
    )

    #plt.bar(
    #    range(len(test_labels)),
    #    np.hstack([test_logit_dict[0], test_logit_dict[1]]),
    #    width=0.7,
    #    color=bar_colors,
    #    align='edge'
    #)

    plt.legend(handles=legend_patches, loc='upper left', ncol=1)

    plt.axhline(y=0.5, xmin=0, xmax=len(test_labels), color='black', linestyle='dashed', linewidth=1.5)

    plt.axvline(x=0.5*len(test_labels), ymin=0, ymax=1, color='black', linestyle='dotted', linewidth=1.5)

    plt.savefig(
        sampler_output_run_paths[i].joinpath('pred_posterior_on_test.png'),
        pil_kwargs={'quality': 100},
        transparent=True,
        bbox_inches='tight',
        pad_inches=0.1
    )

    plt.close()
