# %% Load packages

import numpy as np

# %% Define constants

num_chains = 10

num_epochs = 110000
num_burnin_epochs = 0
diagnostic_iter_thres = 10000

verbose = True
verbose_step = 1000

pred_interval_x_l = -0.5
pred_interval_x_u = 1.5
pred_interval_num = 22
pred_interval_step = (pred_interval_x_u - pred_interval_x_l) / pred_interval_num

pred_interval_x1 = np.linspace(
    pred_interval_x_l + 0.5 * pred_interval_step,
    pred_interval_x_u - 0.5 * pred_interval_step,
    num=pred_interval_num
)
pred_interval_x2 = np.linspace(
    pred_interval_x_l + 0.5 * pred_interval_step,
    pred_interval_x_u - 0.5 * pred_interval_step,
    num=pred_interval_num
)

pred_iter_thres = 100000
