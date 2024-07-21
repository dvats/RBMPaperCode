import torch

from torch.utils.data import DataLoader

from eeyore.datasets import XYDataset

def load_xydataset_from_file(path, dtype=torch.float32, batch_size=None, shuffle=False):
    dataset = XYDataset.from_file(path=path, dtype=dtype)
    dataloader = DataLoader(dataset, batch_size=batch_size or len(dataset), shuffle=shuffle)
    return dataset, dataloader
