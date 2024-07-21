# %% Import packages

import json
import numpy as np

from sklearn.metrics import accuracy_score
from torch.optim import SGD

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.model import model
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.optim.constants import num_epochs
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.optim.dataloaders import training_dataloader
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.optim.sgd.constants import optimizer_output_pilot_path
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.optim.sgd.optimizer import lr, momentum
# from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.optim.sgd.optimizer import loss_fn, lr, momentum
from bnn_mcmc_examples.optim import train

# %% Create output directory if it does not exist

optimizer_output_pilot_path.mkdir(parents=True, exist_ok=True)

# %% Setup SGD optimizer

model.set_params(model.prior.sample())

optimizer = SGD(model.parameters(), lr=lr, momentum=momentum)

# %% Train model

summary = train(
    model,
    training_dataloader,
    optimizer,
    num_epochs,
    # loss_fn=loss_fn,
    save_loss=True,
    save_metric=True,
    save_metric_mean=True,
    pred_fn=lambda labels: labels.squeeze() > 0.5,
    metric_fn=lambda preds, labels: accuracy_score(preds, labels.squeeze())
)

# %% Save solution

with open(optimizer_output_pilot_path.joinpath('solution.csv'), 'w') as file:
    np.savetxt(file, model.get_params().cpu().detach().numpy(), newline='\n', header='')

# %% Save training summary

with open(optimizer_output_pilot_path.joinpath('summary.json'), 'w') as file:
    json.dump(summary, file)
