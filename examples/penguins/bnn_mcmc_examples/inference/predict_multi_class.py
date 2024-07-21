import numpy as np
import torch

def predict_multi_class(num_classes, sample, dataloader, model, pred_prob_path, dropped_samples_path, dtype):
    pred_probs = np.empty([len(dataloader), num_classes])
    nums_dropped_samples = np.empty([len(dataloader), num_classes], dtype=np.int64)

    for i, (x, _) in enumerate(dataloader):
        integral, num_dropped_samples = model.predictive_posterior(sample, x, torch.tensor([[1.]], dtype=dtype))
        pred_probs[i, 1] = integral.item()
        pred_probs[i, 0] = 1. - pred_probs[i, 1]
        nums_dropped_samples[i] = num_dropped_samples

        for j in range(num_classes):
            y = torch.zeros([1, num_classes], dtype=dtype)
            y[0, j] = 1.
            integral, num_dropped_samples = model.predictive_posterior(sample, x, y)
            pred_probs[i, j] = integral.item()
            nums_dropped_samples[i, j] = num_dropped_samples

    np.savetxt(pred_prob_path, pred_probs, delimiter=',')
    np.savetxt(dropped_samples_path, nums_dropped_samples, fmt='%d', delimiter=',')
