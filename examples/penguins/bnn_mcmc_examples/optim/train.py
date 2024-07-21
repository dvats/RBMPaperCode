def train(
    model,
    dataloader,
    optimizer,
    num_epochs,
    loss_fn=None,
    monitor_step=None,
    save_loss=False,
    save_metric=False,
    save_metric_mean=False,
    pred_fn=None,
    metric_fn=None
):
    if loss_fn is None:
        loss_fn = model.loss

    save_summary = save_loss or save_metric or save_metric_mean

    if save_summary and (monitor_step is None):
        monitor_step = len(dataloader)

    loss_vals = []
    metric_vals = []
    metric_mean_vals = []
    num_batches = 0
    num_monitor_batches = 0
    metric_mean_val = 0

    for epoch in range(num_epochs):
        for batch_idx, (input, target) in enumerate(dataloader):
            def closure():
                optimizer.zero_grad()
                output = model(input)
                loss = loss_fn(output, target)
                loss.backward()
                return loss, output

            num_batches = num_batches + 1

            loss_val, output = optimizer.step(closure)

            if save_summary and ((batch_idx + 1) % monitor_step == 0):
                if save_loss:
                    loss_vals.append(loss_val.item())

                if save_metric or save_metric_mean:
                    metric_val = metric_fn(pred_fn(output), target)

                    if save_metric:
                        metric_vals.append(metric_val)

                    if save_metric_mean:
                        num_monitor_batches = num_monitor_batches + 1
                        metric_mean_val = (metric_mean_val * (num_monitor_batches - 1) + metric_val) / num_monitor_batches
                        metric_mean_vals.append(metric_mean_val)

    return {
        'loss_vals': loss_vals,
        'metric_vals': metric_vals,
        'metric_mean_vals': metric_mean_vals,
        'num_batches': num_batches
    }
