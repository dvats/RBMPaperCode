# %% Import packages

from sklearn.metrics import accuracy_score
from torch.optim import SGD

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.model import model
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.optim.constants import (
    num_epochs, num_solutions, verbose, verbose_step
)
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.optim.dataloaders import training_dataloader, test_dataloader
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.optim.sgd.constants import optimizer_output_benchmark_path
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.optim.sgd.optimizer import lr, momentum
# from bnn_mcmc_examples.examples.mlp.noisy_xor.setting2.optim.sgd.optimizer import loss_fn, lr, momentum
from bnn_mcmc_examples.optim import benchmark

# %% Create output directory if it does not exist

optimizer_output_benchmark_path.mkdir(parents=True, exist_ok=True)

# %% Setup SGD optimizer

model.set_params(model.prior.sample())

optimizer = SGD(model.parameters(), lr=lr, momentum=momentum)

# %% Benchmark model

benchmark(
    model,
    training_dataloader,
    optimizer,
    num_solutions,
    num_epochs,
    optimizer_output_benchmark_path,
    # loss_fn=loss_fn,
    validation_loader=test_dataloader,
    pred_fn=lambda labels: labels.squeeze() > 0.5,
    metric_fn=lambda preds, labels: accuracy_score(preds, labels.squeeze()),
    check_fn=lambda acc: acc > 0.85,
    verbose=verbose,
    verbose_step=verbose_step,
    print_runtime=True
)
