import numpy as np
import torch

def predict_binary_class(sample, dataloader, model, pred_prob_path, dropped_samples_path, dtype):
    pred_probs = np.empty([len(dataloader), 2])
    nums_dropped_samples = np.empty([len(dataloader)], dtype=np.int64)

    for i, (x, _) in enumerate(dataloader):
        integral, num_dropped_samples = model.predictive_posterior(sample, x, torch.tensor([[1.]], dtype=dtype))
        pred_probs[i, 1] = integral.item()
        pred_probs[i, 0] = 1. - pred_probs[i, 1]
        nums_dropped_samples[i] = num_dropped_samples

    np.savetxt(pred_prob_path, pred_probs, delimiter=',')
    np.savetxt(dropped_samples_path, nums_dropped_samples, fmt='%d')
