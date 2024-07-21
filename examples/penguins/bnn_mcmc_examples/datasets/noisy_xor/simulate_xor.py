# %% Load packages

import numpy as np

# %% Define function for simulating noisy XOR points

def simulate_xor(n=np.repeat(100, 4), s=0.55):
    x = np.vstack([
        np.column_stack([np.random.rand(n[0]) - s, np.random.rand(n[0]) - s]), # points corresponding to (0, 0)
        np.column_stack([np.random.rand(n[1]) - s, np.random.rand(n[1]) + s]), # points corresponding to (0, 1)
        np.column_stack([np.random.rand(n[2]) + s, np.random.rand(n[2]) - s]), # points corresponding to (1, 0)
        np.column_stack([np.random.rand(n[3]) + s, np.random.rand(n[3]) + s])  # points corresponding to (1, 1)
    ])

    y = np.repeat([0, 1, 1, 0], n)
    
    return x, y
