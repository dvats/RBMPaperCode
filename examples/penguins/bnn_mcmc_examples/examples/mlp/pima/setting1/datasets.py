# %% Import packages

import pandas as pd
import torch

from sklearn.model_selection import train_test_split

from eeyore.constants import torch_to_np_types
from eeyore.datasets import XYDataset

from bnn_mcmc_examples.datasets.pima.data1.constants import data_path
from bnn_mcmc_examples.examples.mlp.pima.setting1.constants import dtype

# %% Load Pima data

x = pd.read_csv(data_path.joinpath('x.csv'))
y = pd.read_csv(data_path.joinpath('y.csv'))

# %% Split data to training and test subsets

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.33, random_state=5000, stratify=y)

# %% Create training dataset

training_dataset = XYDataset(
    torch.from_numpy(x_train.to_numpy(dtype=torch_to_np_types[dtype])),
    torch.from_numpy(y_train.to_numpy(dtype=torch_to_np_types[dtype]))
)

training_dataset.x = \
    (training_dataset.x - torch.mean(training_dataset.x, dim=0, keepdim=True))/ \
    torch.std(training_dataset.x, dim=0, keepdim=True, unbiased=False)

# %% Create test dataset

test_dataset = XYDataset(
    torch.from_numpy(x_test.to_numpy(dtype=torch_to_np_types[dtype])),
    torch.from_numpy(y_test.to_numpy(dtype=torch_to_np_types[dtype]))
)

test_dataset.x = \
    (test_dataset.x - torch.mean(test_dataset.x, dim=0, keepdim=True))/ \
    torch.std(test_dataset.x, dim=0, keepdim=True, unbiased=False)
