
"""
stats.py

One-probe wave statistics using a zero-crossing method.
"""

import numpy as np
from scipy.signal import detrend


def zero_crossing(eta, fs: float, threshold: float | None = None):
    """Zero-crossing analysis of a single wave-probe time series.

    Adapted from classic MATLAB-style implementations (e.g. Urs Neumeier),
    with logic ported to Python/NumPy.

    Parameters
    ----------
    eta : array-like, pandas Series/DataFrame, or 1D NumPy array
        Free-surface elevation [m] as a function of time.
    fs : float
        Sampling frequency [Hz].
    threshold : float, optional
        Minimum crest/trough amplitude to consider a valid wave.
        If None, set to 1% of maximum wave height.

    Returns
    -------
    res : dict
        Dictionary with keys:
        - 'Hs'    : significant wave height (mean of highest 1/3 waves)
        - 'Hmean' : mean wave height
        - 'H1_10' : mean of highest 1/10 waves
        - 'Hmax'  : maximum wave height
        - 'Tmean' : mean wave period
        - 'Ts'    : significant wave period (mean period of highest 1/3 waves)
        - 'wave'  : (N x 2) array [Height, Period] for each detected wave
    names : list[str]
        Names of the primary output parameters.
    """
    import pandas as pd  # local import; only needed if user passes DataFrame/Series

    if fs <= 0:
        raise ValueError("Sampling frequency fs must be greater than zero.")

    # Convert input to 1D NumPy array
    if isinstance(eta, pd.DataFrame):
        data = eta.values
    elif isinstance(eta, pd.Series):
        data = eta.values
    else:
        data = np.asarray(eta)

    if data.ndim > 1:
        if data.shape[1] > 1:
            raise ValueError("zero_crossing: data must be 1D or single-column.")
        data = data[:, 0]

    data = data.astype(float)

    # Invert signal (match original downward-crossing convention)
    data = -data

    # Remove trend/mean
    data = detrend(data)

    nonzero_idx = np.where(data != 0.0)[0]
    names = ['Hs', 'Hmean', 'H1_10', 'Hmax', 'Tmean', 'Ts']

    if len(nonzero_idx) < 2:
        res = {name: np.nan for name in names}
        res['wave'] = np.empty((0, 2))
        return res, names

    d0 = data[nonzero_idx]

    # Sign change → zero crossing
    crossing_mask = d0[:-1] * d0[1:] < 0
    crossing_idx = nonzero_idx[:-1][crossing_mask]

    # Remove first downward crossing if needed
    if data[0] > 0 and len(crossing_idx) > 0:
        crossing_idx = crossing_idx[1:]

    # Keep every second crossing: upward crossings
    crossing_idx = crossing_idx[::2]

    if len(crossing_idx) < 2:
        res = {name: np.nan for name in names}
        res['wave'] = np.empty((0, 2))
        return res, names

    fs = float(fs)
    n_waves = len(crossing_idx) - 1
    wave = np.zeros((n_waves, 4))  # [Height, Crest, Trough, Period]

    for i in range(n_waves):
        start = crossing_idx[i]
        stop = crossing_idx[i + 1]
        segment = data[start:stop]

        crest = np.max(segment)
        trough = -np.min(segment)
        period = (stop - start) / fs

        wave[i, 1] = crest
        wave[i, 2] = trough
        wave[i, 3] = period

    # Total wave height
    wave[:, 0] = wave[:, 1] + wave[:, 2]

    # Threshold
    if threshold is None:
        threshold = 0.01 * np.max(wave[:, 0])
    elif threshold < 0:
        raise ValueError("Wave threshold must not be negative.")

    # Remove/merge small waves
    i = 0
    while i < len(wave):
        small_crest = wave[i, 1] < threshold
        small_trough = wave[i, 2] < threshold

        if small_crest:
            if i > 0:
                wave[i - 1, 1] = max(wave[i - 1, 1], wave[i, 1])
                wave[i - 1, 2] = max(wave[i - 1, 2], wave[i, 2])
                wave[i - 1, 3] += wave[i, 3]
                wave = np.delete(wave, i, axis=0)
                i -= 1
            else:
                wave = np.delete(wave, i, axis=0)
                i -= 1
        elif small_trough:
            if i < len(wave) - 1:
                wave[i, 1] = max(wave[i, 1], wave[i + 1, 1])
                wave[i, 2] = max(wave[i, 2], wave[i + 1, 2])
                wave[i, 3] += wave[i + 1, 3]
                wave = np.delete(wave, i + 1, axis=0)
            else:
                wave = np.delete(wave, i, axis=0)
                i -= 1

        i += 1

    if len(wave) == 0:
        res = {name: np.nan for name in names}
        res['wave'] = np.empty((0, 2))
        return res, names

    wave_unsorted = wave.copy()
    waves_sorted = wave[np.argsort(wave[:, 0])[::-1]]  # sort by height (desc)

    nb = len(waves_sorted)
    n13 = max(1, round(nb / 3))
    n10 = max(1, round(nb / 10))

    Hs = np.mean(waves_sorted[:n13, 0])
    Hmean = np.mean(waves_sorted[:, 0])
    H1_10 = np.mean(waves_sorted[:n10, 0])
    Hmax = np.max(waves_sorted[:, 0])
    Tmean = np.mean(waves_sorted[:, 3])
    Ts = np.mean(waves_sorted[:n13, 3])

    res = {
        'Hs': Hs,
        'Hmean': Hmean,
        'H1_10': H1_10,
        'Hmax': Hmax,
        'Tmean': Tmean,
        'Ts': Ts,
        'wave': wave_unsorted[:, [0, 3]],  # [Height, Period]
    }

    return res, names
