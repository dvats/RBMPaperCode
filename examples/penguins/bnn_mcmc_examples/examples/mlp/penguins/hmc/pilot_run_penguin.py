# %% Import packages

from datetime import timedelta
from timeit import default_timer as timer

from bnn_mcmc_examples.examples.mlp.penguins.constants import num_burnin_epochs, num_epochs, verbose, verbose_step
from bnn_mcmc_examples.examples.mlp.penguins.hmc.constants import sampler_output_pilot_path
from bnn_mcmc_examples.examples.mlp.penguins.hmc.sampler import sampler

# %% Run HMC sampler

start_time = timer()

sampler.run(num_epochs=pow(10,6), num_burnin_epochs=num_burnin_epochs, verbose=verbose, verbose_step=verbose_step)

end_time = timer()
print("Time taken: {}".format(timedelta(seconds=end_time-start_time)))1
# %% Save chain array

sampler.get_chain().to_chainfile(keys=['sample', 'accepted'], path=sampler_output_pilot_path, mode='w')
