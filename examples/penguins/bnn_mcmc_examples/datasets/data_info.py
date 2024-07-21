from pathlib import Path

data_root = Path(__file__).parent.parent.joinpath('data')

data_paths = {name: data_root.joinpath(name) for name in (
    'hawks',
    'MNIST',
    'noisy_xor',
    'penguins',
    'pima'
)}
