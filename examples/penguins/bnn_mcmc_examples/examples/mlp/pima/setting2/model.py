# %% Import packages

import numpy as np
import torch

from torch.distributions import Normal

from eeyore.constants import loss_functions
from eeyore.models import mlp

from bnn_mcmc_examples.examples.mlp.pima.setting2.constants import dtype, mlp_dims, mlp_bias, mlp_activations

# %% Setup MLP model

hparams = mlp.Hyperparameters(dims=mlp_dims, bias=mlp_bias, activations=mlp_activations)

model = mlp.MLP(loss=loss_functions['binary_classification'], hparams=hparams, dtype=dtype)

prior_scale = np.sqrt(10.)

model.prior = Normal(
    torch.zeros(model.num_params(), dtype=model.dtype),
    torch.full([model.num_params()], prior_scale, dtype=model.dtype)
)
