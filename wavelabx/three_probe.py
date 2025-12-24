"""
three_probe.py

Three-gauge array method for incident/reflected separation
following a Goda & Suzuki-style approach.
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend

from .core import compute_wavelength


def _hann_window(N: int) -> np.ndarray:
    """Return Hann window of length N."""
    return np.hanning(N)


def _energy_normalize_window(w: np.ndarray) -> np.ndarray:
    """Normalize window so mean-square equals 1 (preserves variance)."""
    ms = np.mean(w**2)
    return w / np.sqrt(ms) if ms > 0 else w


def three_probe_array(
    eta123: np.ndarray,
    fs: float,
    h: float,
    gpos: tuple[float, float, float],
    plot: bool = False,
    window: str | None = None,
    cond_warn: float = 1e3,
    cond_max: float = 1e6,
    min_retained_energy: float = 0.8,
    figures_dir: str = "figures",
    save_prefix: str = "threeprobe",
    no_bands: int = 5,
) -> dict:
    """Three-gauge array method for separating incident and reflected components
    in random waves.

    Parameters
    ----------
    eta123 : np.ndarray
        (N x 3) array of free-surface elevations [m] from three probes.
        Column 0 = gauge 1, column 1 = gauge 2, column 2 = gauge 3.
    fs : float
        Sampling frequency [Hz].
    h : float
        Water depth [m].
    gpos : tuple[float, float, float]
        Positions [m] of the three gauges along the flume.
    plot : bool, optional
        If True, saves figures (spectra + reconstructed time series) to figures_dir.
    window : str | None, optional
        Window to apply. Supported: None, "hann"/"hanning".
    cond_warn : float
        Print a warning if any pair condition number exceeds this.
    cond_max : float
        Mask frequencies for a pair if condition number exceeds this.
    min_retained_energy : float
        Warn if retained-energy fraction is below this threshold.
    figures_dir : str
        Directory to save figures (created if missing).
    save_prefix : str
        Prefix for saved figure filenames.
    no_bands : int
        Number of bands for band-averaging spectra.

    Returns
    -------
    refanalysis : dict
        Dictionary with (primary keys):
        - 'Hi'    : incident Hm0 estimate
        - 'Hr'    : reflected Hm0 estimate
        - 'Kr'    : reflection coefficient Hr / Hi
        - 'Htot'  : total Hm0 from composite spectrum (gauge 1 spectrum)
        Plus diagnostic fields and spectra arrays.
    """
    eta123 = np.asarray(eta123, dtype=float)
    if eta123.ndim != 2 or eta123.shape[1] != 3:
        raise ValueError("eta123 must be a 2D array with shape (N, 3).")

    # Detrend (demean)
    z = detrend(eta123, axis=0, type="constant")

    # Optional windowing (energy-normalized Hann).
    if isinstance(window, str) and window.strip().lower() in ("none", "off", "false", "0"):
        window = None
    if window is not None:
        if isinstance(window, str) and window.strip().lower() in ("hann", "hanning"):
            w = _energy_normalize_window(_hann_window(z.shape[0]))
        else:
            raise ValueError(f"Unsupported window: {window}")
        z = z * w[:, None]

    N = int(z.shape[0])
    fs = float(fs)
    dt = 1.0 / fs

    # FFT settings
    nfft = N
    df = 1.0 / (nfft * dt)
    half = nfft // 2  # Nyquist excluded by this scheme
    if half < 2:
        raise ValueError("Time series too short for spectral decomposition.")

    # Allocate (for frequencies 1..half-1)
    An = np.zeros((half - 1, 3))
    Bn = np.zeros((half - 1, 3))
    Sn = np.zeros((half - 1, 3))

    # Fourier coefficients and (one-sided) auto-spectrum per gauge
    for j in range(3):
        fn = np.fft.fft(z[:, j], nfft)
        An[:, j] = 2.0 * np.real(fn[1:half]) / nfft
        Bn[:, j] = -2.0 * np.imag(fn[1:half]) / nfft

        fn_sq = fn * np.conj(fn)
        fn_fold = 2.0 * fn_sq[1:half]
        Sn[:, j] = dt * np.real(fn_fold) / nfft

    # Frequency vector for bins 1..half-1
    f = df * np.arange(1, half)

    # Wavenumber from dispersion relation via wavelength
    k = np.zeros_like(f)
    for i, fi in enumerate(f):
        if fi <= 0.0:
            k[i] = 0.0
        else:
            T = 1.0 / fi
            L = compute_wavelength(h, T)
            k[i] = 2.0 * np.pi / L

    # ------------------------------------------------------------
    # Conditioning diagnostic per pair (2x2 complex inversion)
    # ------------------------------------------------------------
    pairs = [(0, 1), (0, 2), (1, 2)]  # (1-2), (1-3), (2-3)
    cond_pair = np.full((len(f), len(pairs)), np.nan, dtype=float)

    for jj, (i1, i2) in enumerate(pairs):
        x1, x2 = gpos[i1], gpos[i2]
        for ii in range(len(f)):
            if f[ii] <= 0.0 or k[ii] <= 0.0:
                continue
            M = np.array(
                [
                    [np.exp(-1j * k[ii] * x1), np.exp(1j * k[ii] * x1)],
                    [np.exp(-1j * k[ii] * x2), np.exp(1j * k[ii] * x2)],
                ],
                dtype=complex,
            )
            try:
                cond_pair[ii, jj] = np.linalg.cond(M)
            except Exception:
                cond_pair[ii, jj] = np.inf

    if np.any(cond_pair > cond_warn):
        worst = np.nanmax(cond_pair)
        print(
            f"[WaveLabX] Warning: three-probe pair inversions are ill-conditioned at some "
            f"frequencies (max cond ≈ {worst:.2e})."
        )

    bad_cond_pair = cond_pair > cond_max

    # Pair indexing mapping to columns 0,1,2 -> (1-2),(1-3),(2-3)
    g1 = [0, 0, 1]
    g2 = [1, 2, 2]

    Ainc = np.zeros((len(k), 3))
    Binc = np.zeros((len(k), 3))
    Aref = np.zeros((len(k), 3))
    Bref = np.zeros((len(k), 3))

    nmin = np.zeros(3, dtype=int)
    nmax = np.zeros(3, dtype=int)

    # Compute per-pair incident/reflected Fourier coefficients
    for j in range(3):
        A1 = An[:, g1[j]]
        A2 = An[:, g2[j]]
        B1 = Bn[:, g1[j]]
        B2 = Bn[:, g2[j]]

        pos1 = gpos[g1[j]]
        pos2 = gpos[g2[j]]
        dx = abs(pos2 - pos1)

        term1 = -A2 * np.sin(k * pos1) + A1 * np.sin(k * pos2) + B2 * np.cos(k * pos1) - B1 * np.cos(k * pos2)
        term2 =  A2 * np.cos(k * pos1) - A1 * np.cos(k * pos2) + B2 * np.sin(k * pos1) - B1 * np.sin(k * pos2)
        term3 = -A2 * np.sin(k * pos1) + A1 * np.sin(k * pos2) - B2 * np.cos(k * pos1) + B1 * np.cos(k * pos2)
        term4 =  A2 * np.cos(k * pos1) - A1 * np.cos(k * pos2) - B2 * np.sin(k * pos1) + B1 * np.sin(k * pos2)

        denom = 2.0 * np.sin(k * dx + 1e-16)

        Ainc[:, j] = term1 / denom
        Binc[:, j] = term2 / denom
        Aref[:, j] = term3 / denom
        Bref[:, j] = term4 / denom

        # Mask ill-conditioned frequencies for this pair
        if np.any(bad_cond_pair[:, j]):
            Ainc[bad_cond_pair[:, j], j] = np.nan
            Binc[bad_cond_pair[:, j], j] = np.nan
            Aref[bad_cond_pair[:, j], j] = np.nan
            Bref[bad_cond_pair[:, j], j] = np.nan

        # Valid wavelength band (Goda guideline) for this pair
        Lmin = dx / 0.45
        Lmax = dx / 0.05
        kmax = 2.0 * np.pi / Lmin
        kmin = 2.0 * np.pi / Lmax
        ind = np.where((k <= kmax) & (k >= kmin))[0]

        if len(ind) == 0:
            nmin[j] = 0
            nmax[j] = len(k) - 1
        else:
            nmin[j] = int(ind[0])
            nmax[j] = int(ind[-1])

    # Mask outside valid k-range per pair
    idx = np.arange(len(k))
    for j in range(3):
        mask = (idx < nmin[j]) | (idx > nmax[j])
        Ainc[mask, j] = np.nan
        Binc[mask, j] = np.nan
        Aref[mask, j] = np.nan
        Bref[mask, j] = np.nan

    # Global usable range (union of pair-valid ranges)
    global_min = int(np.min(nmin))
    global_max = int(np.max(nmax))
    rng = np.arange(global_min, global_max + 1, dtype=int)

    # Retained energy fraction diagnostic (based on gauge-1 spectrum)
    S_total = Sn[:, 0].copy()
    if len(S_total) > 0:
        S_total[0] = 0.0
    total_energy = np.nansum(S_total) * df
    retained_energy = np.nansum(S_total[rng]) * df if total_energy > 0 else 0.0
    retained_energy_fraction = (retained_energy / total_energy) if total_energy > 0 else np.nan

    if (not np.isnan(retained_energy_fraction)) and (retained_energy_fraction < min_retained_energy):
        print(
            f"[WaveLabX] Warning: three-probe method retained only {retained_energy_fraction:.2%} "
            f"of spectral energy (threshold {min_retained_energy:.0%}). Results may be unreliable."
        )

    # Average across valid pairs at each frequency (nanmean handles masked values)
    Ainc_av = np.zeros_like(k)
    Binc_av = np.zeros_like(k)
    Aref_av = np.zeros_like(k)
    Bref_av = np.zeros_like(k)

    Ainc_av[rng] = np.nanmean(Ainc[rng, :], axis=1)
    Binc_av[rng] = np.nanmean(Binc[rng, :], axis=1)
    Aref_av[rng] = np.nanmean(Aref[rng, :], axis=1)
    Bref_av[rng] = np.nanmean(Bref[rng, :], axis=1)

    # Spectra from averaged Fourier coefficients
    Si = (Ainc_av[rng] ** 2 + Binc_av[rng] ** 2) / (2.0 * df)
    Sr = (Aref_av[rng] ** 2 + Bref_av[rng] ** 2) / (2.0 * df)
    Ssum = Si + Sr  # what we plot as "composite"
    Sf = Sn[rng, 0]  # gauge 1 spectrum, used for Htot only
    flim = f[rng]

    Ei = np.nansum(Si) * df
    Er = np.nansum(Sr) * df
    mo = np.nansum(Sf) * df

    Htot = 4.0 * np.sqrt(mo)
    Hi = 4.0 * np.sqrt(Ei)
    Hr = 4.0 * np.sqrt(Er)
    refco = Hr / (Hi + 1e-16)

    # ----------------------------------------------------------------------
    # Plotting (NO measured/gauge spectrum in plots; save to figures/)
    # ----------------------------------------------------------------------
    if plot:
        os.makedirs(figures_dir, exist_ok=True)

        plt.rcParams.update(
            {
                "font.size": 9,
                "axes.labelsize": 9,
                "axes.titlesize": 10,
                "legend.fontsize": 8,
                "xtick.labelsize": 8,
                "ytick.labelsize": 8,
                "figure.dpi": 150,
                "savefig.dpi": 300,
            }
        )

        # Band-averaging (simple)
        no_bands = max(1, int(no_bands))
        n_pts = len(Si)
        band_len = max(1, n_pts // no_bands)

        flim_band, Si_band, Sr_band, Ssum_band = [], [], [], []
        for j in range(0, n_pts, band_len):
            j_end = min(j + band_len, n_pts)
            flim_band.append(np.mean(flim[j:j_end]))
            Si_band.append(np.mean(Si[j:j_end]))
            Sr_band.append(np.mean(Sr[j:j_end]))
            Ssum_band.append(np.mean(Ssum[j:j_end]))

        flim_band = np.asarray(flim_band)
        Si_band = np.asarray(Si_band)
        Sr_band = np.asarray(Sr_band)
        Ssum_band = np.asarray(Ssum_band)

        # (1) Spectra figure (3 panels)
        fig, axes = plt.subplots(3, 1, figsize=(7, 8), constrained_layout=True)

        ax = axes[0]
        ax.plot(flim, Si, label="Incident", linewidth=1.2)
        ax.plot(flim, Sr, label="Reflected", linestyle="--", linewidth=1.2)
        ax.plot(flim, Ssum, label="Incident + Reflected", linewidth=1.0)
        ax.set_title("Incident, Reflected, and Composite Spectra")
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(r"$S(f)$ [m$^2$s]")
        ax.grid(alpha=0.3)
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax = axes[1]
        ax.plot(flim_band, Si_band, label="Incident (band avg)", linestyle=":")
        ax.plot(flim_band, Sr_band, label="Reflected (band avg)", linestyle="-.")
        ax.plot(flim_band, Ssum_band, label="Incident + Reflected (band avg)")
        ax.set_title("Band-Averaged Spectra")
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(r"$S(f)$ [m$^2$s]")
        ax.grid(alpha=0.3)
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax = axes[2]
        ax.plot(flim_band, Si_band, label="Incident (band avg)", linewidth=1.2)
        ax.set_title("Incident Band-Averaged Spectrum")
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(r"$S(f)$ [m$^2$s]")
        ax.grid(alpha=0.3)
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.savefig(os.path.join(figures_dir, f"{save_prefix}_spectra.png"), bbox_inches="tight")

        # (2) Reconstructed incident/reflected time series (optional)
        t = np.arange(N) * dt
        eta_inc = np.zeros(N)
        eta_ref = np.zeros(N)

        # Only sum valid frequencies (rng); outside rng are zeros by construction
        for i in rng:
            omega = 2.0 * np.pi * f[i]
            eta_inc += Ainc_av[i] * np.cos(omega * t) + Binc_av[i] * np.sin(omega * t)
            eta_ref += Aref_av[i] * np.cos(omega * t) + Bref_av[i] * np.sin(omega * t)

        fig2, axes2 = plt.subplots(2, 1, figsize=(7, 5), constrained_layout=True)

        ax = axes2[0]
        ax.plot(t, eta_inc, linewidth=1.0)
        ax.set_title("Reconstructed Incident Wave Time Series")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(r"$\eta_{\mathrm{inc}}$ [m]")
        ax.grid(alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax = axes2[1]
        ax.plot(t, eta_ref, linewidth=1.0)
        ax.set_title("Reconstructed Reflected Wave Time Series")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(r"$\eta_{\mathrm{ref}}$ [m]")
        ax.grid(alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig2.savefig(os.path.join(figures_dir, f"{save_prefix}_incident_reflected.png"), bbox_inches="tight")

    return {
        # Primary output keys
        "Kr": refco,
        "Hi": Hi,
        "Hr": Hr,
        "Htot": Htot,
        # Aliases
        "refco": refco,
        "Htotal": Htot,
        # Spectra and frequencies
        "Si": Si,
        "Sr": Sr,
        "Ssum": Ssum,
        "Sf": Sf,  # (gauge spectrum; kept for diagnostics)
        "f": flim,
        # Diagnostics
        "cond_pair": cond_pair,
        "bad_cond_pair": bad_cond_pair,
        "retained_energy_fraction": retained_energy_fraction,
        "valid_index_range": (int(global_min), int(global_max)),
    }
