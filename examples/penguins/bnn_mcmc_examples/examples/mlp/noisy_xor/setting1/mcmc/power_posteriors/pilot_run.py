# %% Import packages

from datetime import timedelta
from timeit import default_timer as timer

from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.constants import (
    num_burnin_epochs, num_epochs, verbose, verbose_step
)
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.power_posteriors.constants import sampler_output_pilot_path
from bnn_mcmc_examples.examples.mlp.noisy_xor.setting1.mcmc.power_posteriors.sampler import sampler

# %% Run power posterior sampler

start_time = timer()

sampler.run(num_epochs=num_epochs, num_burnin_epochs=num_burnin_epochs, verbose=verbose, verbose_step=verbose_step)

end_time = timer()
print("Time taken: {}".format(timedelta(seconds=end_time-start_time)))

# %% Save chain array

sampler.get_chain().to_chainfile(keys=['sample'], path=sampler_output_pilot_path, mode='w')
