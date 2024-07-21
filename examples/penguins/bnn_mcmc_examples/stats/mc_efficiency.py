def mc_efficiency(
    paths,
    keys=['rhat', 'ess', 'runtime', 'efficiency', 'rel_efficiency'],
    rhat_fname='multi_rhat.txt',
    ess_fname='mean_multi_ess.txt',
    runtime_fname='mean_runtime.txt'
):
    num_paths = len(paths)
    compute_efficiency = any(k in keys for k in ['efficiency', 'rel_efficiency'])
    all_keys = keys.copy()
    if compute_efficiency:
        for key in ['ess', 'runtime']:
            if key not in keys:
                all_keys.append(key)

    summary = {key: [] for key in all_keys}
    fnames = {'rhat': rhat_fname, 'ess': ess_fname, 'runtime': runtime_fname}
    parsers = {'rhat': lambda x: float(x), 'ess': lambda x: int(float(x)), 'runtime': lambda x: float(x)}

    for key in [read_key for read_key in all_keys if read_key not in ['efficiency', 'rel_efficiency']]:
        for path in paths:
            with open(path.joinpath(fnames[key]), 'r') as file:
                summary[key].append(parsers[key](file.readline().rstrip()))

    if compute_efficiency:        
        for i in range(num_paths):
            summary['efficiency'].append(summary['ess'][i] / summary['runtime'][i])

        if 'rel_efficiency' in keys:
            summary['rel_efficiency'].append(1.)
        
            for i in range(1, num_paths):
                summary['rel_efficiency'].append(summary['efficiency'][i] / summary['efficiency'][0])

    return {key: value for key, value in summary.items() if key in keys}
