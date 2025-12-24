
"""
three_probe.py

Three-gauge array method for incident/reflected separation
following a Goda & Suzuki-style approach.
"""

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
        If True, plots spectra and reconstructed time series.

    Returns
    -------
    refanalysis : dict
        Dictionary with:
        - 'Hi'    : incident Hm0 estimate
        - 'Hr'    : reflected Hm0 estimate
        - 'refco' : reflection coefficient Hr / Hi
        - 'Htot'  : total Hm0 from composite spectrum
    """
    eta123 = np.asarray(eta123, dtype=float)
    if eta123.ndim != 2 or eta123.shape[1] != 3:
        raise ValueError("eta123 must be a 2D array with shape (N, 3).")

    z = detrend(eta123, axis=0, type="constant")

    # Optional windowing (energy-normalized Hann).
    # Accept common ways users may disable it: None, "none", "off", "false".
    if isinstance(window, str) and window.strip().lower() in ("none", "off", "false", "0"):
        window = None
    if window is not None:
        if isinstance(window, str) and window.lower() in ("hann", "hanning"):
            w = _energy_normalize_window(_hann_window(z.shape[0]))
        else:
            raise ValueError(f"Unsupported window: {window}")
        z = z * w[:, None]
    N = z.shape[0]
    fs = float(fs)
    dt = 1.0 / fs

    # FFT settings
    nfft = N
    df = 1.0 / (nfft * dt)
    half = nfft // 2  # excluding Nyquist in this scheme

    An = np.zeros((half - 1, 3))
    Bn = np.zeros((half - 1, 3))
    Sn = np.zeros((half - 1, 3))

    # Fourier coefficients for each gauge
    for j in range(3):
        fn = np.fft.fft(z[:, j], nfft)
        An[:, j] = 2.0 * np.real(fn[1:half]) / nfft
        Bn[:, j] = -2.0 * np.imag(fn[1:half]) / nfft
        fn_sq = fn * np.conj(fn)
        fn_fold = 2.0 * fn_sq[1:half]
        Sn[:, j] = dt * np.real(fn_fold) / nfft

    # Frequency vector for 1..half-1
    f = df * np.arange(1, half)

    # Wavenumber from dispersion
    k = np.zeros_like(f)
    for i, fi in enumerate(f):
        if fi <= 0.0:
            k[i] = 0.0
        else:
            T = 1.0 / fi
            L = compute_wavelength(h, T)
            k[i] = 2.0 * np.pi / L

    # Gauge pairs (0-based)

    # ------------------------------------------------------------
    # Conditioning diagnostic per pair (2x2 complex inversion)
    # ------------------------------------------------------------
    pairs = [(0, 1), (0, 2), (1, 2)]  # corresponds to (1-2), (1-3), (2-3)
    cond_pair = np.full((len(f), 3), np.nan, dtype=float)
    for jj, (i1, i2) in enumerate(pairs):
        x1, x2 = gpos[i1], gpos[i2]
        for ii in range(len(f)):
            if f[ii] <= 0.0 or k[ii] <= 0.0:
                continue
            M = np.array([
                [np.exp(-1j * k[ii] * x1), np.exp( 1j * k[ii] * x1)],
                [np.exp(-1j * k[ii] * x2), np.exp( 1j * k[ii] * x2)],
            ], dtype=complex)
            try:
                cond_pair[ii, jj] = np.linalg.cond(M)
            except Exception:
                cond_pair[ii, jj] = np.inf

    if np.any(cond_pair > cond_warn):
        worst = np.nanmax(cond_pair)
        print(f"[WaveLabX] Warning: three-probe pair inversions are ill-conditioned at some frequencies (max cond ≈ {worst:.2e}).")
    bad_cond_pair = cond_pair > cond_max
    g1 = [0, 0, 1]
    g2 = [1, 2, 2]

    Ainc = np.zeros_like(An)
    Binc = np.zeros_like(An)
    Aref = np.zeros_like(An)
    Bref = np.zeros_like(An)

    nmin = np.zeros(3, dtype=int)
    nmax = np.zeros(3, dtype=int)

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


        # Valid wavelength range for this pair
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
    for j in range(3):
        mask = (np.arange(len(k)) < nmin[j]) | (np.arange(len(k)) > nmax[j])
        Ainc[mask, j] = np.nan
        Binc[mask, j] = np.nan
        Aref[mask, j] = np.nan
        Bref[mask, j] = np.nan

    global_min = int(np.min(nmin))
    global_max = int(np.max(nmax))
    rng = np.arange(global_min, global_max + 1, dtype=int)

    # Retained energy fraction diagnostic (based on composite spectrum of gauge 1)
    S_total = Sn[:, 0].copy()
    if len(S_total) > 0:
        S_total[0] = 0.0
    total_energy = np.nansum(S_total) * df
    retained_energy = np.nansum(S_total[rng]) * df if total_energy > 0 else 0.0
    retained_energy_fraction = (retained_energy / total_energy) if total_energy > 0 else np.nan
    if (not np.isnan(retained_energy_fraction)) and (retained_energy_fraction < min_retained_energy):
        print(f"[WaveLabX] Warning: three-probe method retained only {retained_energy_fraction:.2%} of spectral energy (threshold {min_retained_energy:.0%}). Results may be unreliable.")

    # (The diagnostic above already computes retained_energy_fraction and warns.)

    Ainc_av = np.zeros_like(k)
    Binc_av = np.zeros_like(k)
    Aref_av = np.zeros_like(k)
    Bref_av = np.zeros_like(k)

    Ainc_av[rng] = np.nanmean(Ainc[rng, :], axis=1)
    Binc_av[rng] = np.nanmean(Binc[rng, :], axis=1)
    Aref_av[rng] = np.nanmean(Aref[rng, :], axis=1)
    Bref_av[rng] = np.nanmean(Bref[rng, :], axis=1)

    # Spectra
    Si = (Ainc_av[rng]**2 + Binc_av[rng]**2) / (2.0 * df)
    Sr = (Aref_av[rng]**2 + Bref_av[rng]**2) / (2.0 * df)
    Sf = Sn[rng, 0]  # total spectra from gauge 1
    flim = f[rng]

    Ei = np.sum(Si) * df
    Er = np.sum(Sr) * df
    mo = np.sum(Sf) * df

    Htot = 4.0 * np.sqrt(mo)
    Hi = 4.0 * np.sqrt(Ei)
    Hr = 4.0 * np.sqrt(Er)
    refco = Hr / (Hi + 1e-16)

    # if plot:
    #     # Band-averaged spectra
    #     no_bands = 5
    #     n_pts = len(Si)
    #     band_len = max(1, n_pts // no_bands)

    #     flim = f[rng]
    #     flim_band = []
    #     Si_band = []
    #     Sr_band = []
    #     Sf_band = []

    #     for j in range(0, n_pts, band_len):
    #         j_end = min(j + band_len, n_pts)
    #         flim_band.append(np.mean(flim[j:j_end]))
    #         Si_band.append(np.mean(Si[j:j_end]))
    #         Sr_band.append(np.mean(Sr[j:j_end]))
    #         Sf_band.append(np.mean(Sf[j:j_end]))

    #     flim_band = np.array(flim_band)
    #     Si_band = np.array(Si_band)
    #     Sr_band = np.array(Sr_band)
    #     Sf_band = np.array(Sf_band)

    #     plt.figure(figsize=(8, 10))

    #     # Detailed spectra
    #     plt.subplot(3, 1, 1)
    #     plt.plot(flim, Si, ':b', label="Incident")
    #     plt.plot(flim, Sr, 'r-.', label="Reflected")
    #     plt.plot(flim, Sf, 'k', label="Composite")
    #     plt.xlabel("Frequency [Hz]")
    #     plt.ylabel("S(f) [m²·s]")
    #     plt.title("Incident, Reflected, and Composite Spectra")
    #     plt.legend()
    #     plt.grid(True)

    #     # Band-averaged spectra
    #     plt.subplot(3, 1, 2)
    #     plt.plot(flim_band, Si_band, ':b', label="Incident (band avg)")
    #     plt.plot(flim_band, Sr_band, 'r-.', label="Reflected (band avg)")
    #     plt.plot(flim_band, Sf_band, 'k', label="Composite (band avg)")
    #     plt.xlabel("Frequency [Hz]")
    #     plt.ylabel("S(f) [m²·s]")
    #     plt.title("Band-Averaged Spectra")
    #     plt.legend()
    #     plt.grid(True)

    #     # Incident band-averaged spectrum
    #     plt.subplot(3, 1, 3)
    #     plt.plot(flim_band, Si_band, '-b')
    #     plt.xlabel("Frequency [Hz]")
    #     plt.ylabel("S(f) [m²·s]")
    #     plt.title("Incident Band-Averaged Spectrum")
    #     plt.grid(True)

    #     plt.tight_layout()

    #     # Reconstruct incident/reflected time series
    #     t = np.arange(N) * dt
    #     eta_inc = np.zeros(N)
    #     eta_ref = np.zeros(N)

    #     for i in range(len(k)):
    #         omega = 2.0 * np.pi * f[i]
    #         eta_inc += Ainc_av[i] * np.cos(omega * t) + Binc_av[i] * np.sin(omega * t)
    #         eta_ref += Aref_av[i] * np.cos(omega * t) + Bref_av[i] * np.sin(omega * t)

    #     plt.figure(figsize=(8, 6))
    #     plt.subplot(2, 1, 1)
    #     plt.plot(t, eta_inc, 'b')
    #     plt.xlabel("Time [s]")
    #     plt.ylabel("η_inc [m]")
    #     plt.title("Reconstructed Incident Wave Time Series")
    #     plt.grid(True)

    #     plt.subplot(2, 1, 2)
    #     plt.plot(t, eta_ref, 'r')
    #     plt.xlabel("Time [s]")
    #     plt.ylabel("η_ref [m]")
    #     plt.title("Reconstructed Reflected Wave Time Series")
    #     plt.grid(True)

    #     plt.tight_layout()


        # ----------------------------------------------------------------------
    # Enhanced, publication-quality plotting
    # ----------------------------------------------------------------------
    if plot:
        plt.rcParams.update({
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.dpi": 150,
            "savefig.dpi": 300,    # high resolution for journals
        })

        # -------------------------------
        # Prepare band-averaged spectra
        # -------------------------------
        no_bands = 5
        n_pts = len(Si)
        band_len = max(1, n_pts // no_bands)

        flim = f[rng]
        flim_band = []
        Si_band = []
        Sr_band = []
        Sf_band = []

        for j in range(0, n_pts, band_len):
            j_end = min(j + band_len, n_pts)
            flim_band.append(np.mean(flim[j:j_end]))
            Si_band.append(np.mean(Si[j:j_end]))
            Sr_band.append(np.mean(Sr[j:j_end]))
            Sf_band.append(np.mean(Sf[j:j_end]))

        flim_band = np.array(flim_band)
        Si_band  = np.array(Si_band)
        Sr_band  = np.array(Sr_band)
        Sf_band  = np.array(Sf_band)

        # -------------------------------
        # (1) Detailed + Band-Averaged Spectra
        # -------------------------------
        fig1, axes = plt.subplots(3, 1, figsize=(7, 8), constrained_layout=True)

        # --- Panel 1: detailed spectra
        ax = axes[0]
        ax.plot(flim, Si, label="Incident", color="#1f77b4", linewidth=1.2)
        ax.plot(flim, Sr, label="Reflected", color="#d62728", linestyle="--", linewidth=1.2)
        ax.plot(flim, Sf, label="Composite", color="black", linewidth=1.0)

        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(r"$S(f)$ [m²·s]")
        ax.set_title("Incident, Reflected, and Composite Spectra")
        ax.grid(alpha=0.3)
        ax.legend()

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        # --- Panel 2: band-averaged
        ax = axes[1]
        ax.plot(flim_band, Si_band, label="Incident (band avg)", color="#1f77b4", linestyle=":")
        ax.plot(flim_band, Sr_band, label="Reflected (band avg)", color="#d62728", linestyle="-.")
        ax.plot(flim_band, Sf_band, label="Composite (band avg)", color="black")

        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(r"$S(f)$ [m²·s]")
        ax.set_title("Band-Averaged Spectra")
        ax.grid(alpha=0.3)
        ax.legend()

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        # --- Panel 3: incident band-averaged
        ax = axes[2]
        ax.plot(flim_band, Si_band, color="#1f77b4", linewidth=1.2)
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(r"$S(f)$ [m²·s]")
        ax.set_title("Incident Band-Averaged Spectrum")
        ax.grid(alpha=0.3)

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        # Save detailed spectrum figure
        fig1.savefig("threeprobe_spectra.png", dpi=300, bbox_inches="tight")
        # fig1.savefig("threeprobe_spectra.pdf", bbox_inches="tight")

        # -------------------------------
        # (2) Reconstructed Incident / Reflected Time Series
        # -------------------------------
        t = np.arange(N) * dt
        eta_inc = np.zeros(N)
        eta_ref = np.zeros(N)

        for i in range(len(k)):
            omega = 2.0 * np.pi * f[i]
            eta_inc += Ainc_av[i] * np.cos(omega * t) + Binc_av[i] * np.sin(omega * t)
            eta_ref += Aref_av[i] * np.cos(omega * t) + Bref_av[i] * np.sin(omega * t)

        fig2, axes2 = plt.subplots(2, 1, figsize=(7, 5), constrained_layout=True)

        # Incident
        ax = axes2[0]
        ax.plot(t, eta_inc, color="#1f77b4", linewidth=1.0)
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(r"$\eta_{\mathrm{inc}}$ [m]")
        ax.set_title("Reconstructed Incident Wave Time Series")
        ax.grid(alpha=0.3)

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        # Reflected
        ax = axes2[1]
        ax.plot(t, eta_ref, color="#d62728", linewidth=1.0)
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(r"$\eta_{\mathrm{ref}}$ [m]")
        ax.set_title("Reconstructed Reflected Wave Time Series")
        ax.grid(alpha=0.3)

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        fig2.savefig("incident_reflected.png", dpi=300, bbox_inches="tight")
        # fig2.savefig("threeprobe_timeseries.pdf", bbox_inches="tight")


    return {
        # Primary output keys
        'Kr': refco,
        'Hi': Hi,
        'Hr': Hr,
        'Htot': Htot,

        # Backward/reader-friendly aliases
        'refco': refco,
        'Htotal': Htot,
        'Si': Si,
        'Sr': Sr,
        'Sf': Sf,
        'f': flim,
        'cond_pair': cond_pair,
        'bad_cond_pair': bad_cond_pair,
        'retained_energy_fraction': retained_energy_fraction,
        'valid_index_range': (int(global_min), int(global_max)),
    }