"""
two_probe.py

Two-probe reflection analysis using a Goda/Suzuki-style spectral method,
with scaling consistent with the three-probe array implementation and
wavenumbers computed via core.compute_wavelength.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend

from .core import GRAVITY, compute_wavelength

def _hann_window(N: int) -> np.ndarray:
    """Return Hann window of length N."""
    return np.hanning(N)

def _energy_normalize_window(w: np.ndarray) -> np.ndarray:
    """Normalize window so mean-square equals 1 (preserves variance)."""
    ms = np.mean(w**2)
    return w / np.sqrt(ms) if ms > 0 else w

def two_probe_goda(
    eta12: np.ndarray,
    fs: float,
    h: float,
    gpos: tuple[float, float],
    plot: bool = False,
    window: str | None = None,
    cond_warn: float = 1e3,
    cond_max: float = 1e6,
) -> dict:
    """Two-wave-probe reflection analysis (Goda-style method).

    Parameters
    ----------
    eta12 : np.ndarray
        (N x 2) array of free-surface elevations [m] from two probes:
        column 0 = probe 1, column 1 = probe 2.
    fs : float
        Sampling frequency [Hz].
    h : float
        Mean water depth [m].
    gpos : tuple of float
        (x1, x2) positions of the two probes [m] along the flume.
    plot : bool, optional
        If True, plots spectral incident/reflected Hm0 vs frequency.

    Returns
    -------
    refanalysis : dict
        Dictionary containing:
        - 'Kr' : reflection coefficient Kr = Hr / Hi
        - 'Hi' : incident Hm0 [m]
        - 'Hr' : reflected Hm0 [m]
        - 'Si' : incident spectral density array [m^2/Hz]
        - 'Sr' : reflected spectral density array [m^2/Hz]
        - 'f'  : frequency array [Hz]
    """
    # ------------------------------------------------------------
    # Input checks and detrending
    # ------------------------------------------------------------
    eta12 = np.asarray(eta12, dtype=float)
    if eta12.ndim != 2 or eta12.shape[1] != 2:
        raise ValueError("eta12 must be a 2D array with shape (N, 2).")

    eta1 = detrend(eta12[:, 0], type="constant")
    eta2 = detrend(eta12[:, 1], type="constant")

    N = min(len(eta1), len(eta2))
    eta1 = eta1[:N]
    eta2 = eta2[:N]

    # Optional windowing (energy-normalized Hann).
    # Accept common ways users may disable it: None, "none", "off", "false".
    if isinstance(window, str) and window.strip().lower() in ("none", "off", "false", "0"):
        window = None
    if window is not None:
        if isinstance(window, str) and window.lower() in ("hann", "hanning"):
            w = _energy_normalize_window(_hann_window(N))
        else:
            raise ValueError(f"Unsupported window: {window}")
        eta1 = eta1 * w
        eta2 = eta2 * w

    fs = float(fs)
    dt = 1.0 / fs
    Nfft = N

    # ------------------------------------------------------------
    # Frequency vector and FFT (real-valued rfft)
    # ------------------------------------------------------------
    f = np.fft.rfftfreq(Nfft, d=dt)  # 0 .. fs/2
    if f.size > 1:
        df = f[1] - f[0]
    else:
        df = fs / Nfft

    # FFT with (2/N) scaling (same spirit as MATLAB convention)
    A1 = np.fft.rfft(eta1, n=Nfft) * (2.0 / Nfft)
    A2 = np.fft.rfft(eta2, n=Nfft) * (2.0 / Nfft)

    Re1, Im1 = A1.real, A1.imag
    Re2, Im2 = A2.real, A2.imag

    # ------------------------------------------------------------
    # Wavenumber k(f) from dispersion using compute_wavelength
    # ------------------------------------------------------------
    k = np.zeros_like(f)
    for i in range(1, len(f)):  # skip f=0
        T = 1.0 / f[i]
        L = compute_wavelength(h, T)
        if L > 0.0:
            k[i] = 2.0 * np.pi / L
        else:
            k[i] = 0.0
    # k[0] stays 0.0

    dx = abs(gpos[1] - gpos[0])

    # ------------------------------------------------------------
    # Conditioning diagnostic (2x2 complex inversion matrix)
    # ------------------------------------------------------------
    x1, x2 = gpos
    cond_M = np.full_like(f, np.nan, dtype=float)
    for ii in range(len(f)):
        if f[ii] <= 0.0 or k[ii] <= 0.0:
            continue
        M = np.array([
            [np.exp(-1j * k[ii] * x1), np.exp( 1j * k[ii] * x1)],
            [np.exp(-1j * k[ii] * x2), np.exp( 1j * k[ii] * x2)],
        ], dtype=complex)
        try:
            cond_M[ii] = np.linalg.cond(M)
        except Exception:
            cond_M[ii] = np.inf

    # Mask frequencies with extreme ill-conditioning
    bad_cond = cond_M > cond_max
    if np.any(cond_M > cond_warn):
        # lightweight warning (no logging dependency)
        worst = np.nanmax(cond_M)
        print(f"[WaveLabX] Warning: two-probe inversion is ill-conditioned at some frequencies (max cond ≈ {worst:.2e}).")

    # Amplification factor C = 1 / (2 |sin(k dx)|)
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = 2.0 * np.abs(np.sin(k * dx))
        C = 1.0 / denom
    C[np.isnan(C) | np.isinf(C)] = 1.0
    if C.size > 0:
        C[0] = 1.0
    C[C > 1.0] = 1.0

    cos_kdx = np.cos(k * dx)
    sin_kdx = np.sin(k * dx)

    # ------------------------------------------------------------
    # Incident and reflected amplitudes (Goda algebra)
    # ------------------------------------------------------------
    a_inc = C * np.sqrt(
        (Re2 - Re1 * cos_kdx + Im1 * sin_kdx) ** 2
        + (Im2 - Re1 * sin_kdx - Im1 * cos_kdx) ** 2
    )

    a_ref = C * np.sqrt(
        (Re2 - Re1 * cos_kdx - Im1 * sin_kdx) ** 2
        + (Im2 + Re1 * sin_kdx - Im1 * cos_kdx) ** 2
    )

    # ------------------------------------------------------------
    # Spectral densities (consistent with three-probe method)
    # ------------------------------------------------------------
    # S(f) = a^2 / (2 df)  [m^2/Hz]
    Si = (a_inc ** 2) / (2.0 * df)
    Sr = (a_ref ** 2) / (2.0 * df)

    # Keep a copy for diagnostics (before masking)
    Si_raw = Si.copy()
    Sr_raw = Sr.copy()

    # Apply conditioning mask
    if np.any(bad_cond):
        Si = Si.copy(); Sr = Sr.copy()
        Si[bad_cond] = np.nan
        Sr[bad_cond] = np.nan

    # Remove DC from energy integration
    Si_int = Si.copy()
    Sr_int = Sr.copy()
    if Si_int.size > 0:
        Si_int[0] = 0.0
        Sr_int[0] = 0.0

    Ei = np.nansum(Si_int) * df  # incident variance-like energy
    Er = np.nansum(Sr_int) * df  # reflected variance-like energy

    # Retained-energy diagnostic (how much usable energy remains after masking)
    Si_raw_int = Si_raw.copy()
    Sr_raw_int = Sr_raw.copy()
    if Si_raw_int.size > 0:
        Si_raw_int[0] = 0.0
        Sr_raw_int[0] = 0.0
    Ei_raw = np.nansum(Si_raw_int) * df
    Er_raw = np.nansum(Sr_raw_int) * df
    retained_energy_fraction_i = (Ei / Ei_raw) if Ei_raw > 0 else np.nan
    retained_energy_fraction_r = (Er / Er_raw) if Er_raw > 0 else np.nan

    Hi = 4.004 * np.sqrt(Ei) if Ei > 0 else 0.0
    Hr = 4.004 * np.sqrt(Er) if Er > 0 else 0.0
    Kr = Hr / Hi if Hi > 0 else np.nan

    # ------------------------------------------------------------
    # Optional plotting
    # ------------------------------------------------------------
    # if plot:
    #     plt.figure()
    #     plt.plot(f, 4.0 * np.sqrt(Si), label="Incident Hm0(f)")
    #     plt.plot(f, 4.0 * np.sqrt(Sr), "--", label="Reflected Hm0(f)")
    #     plt.xlabel("Frequency [Hz]")
    #     plt.ylabel("S(f) [m²·s]")
    #     plt.title("Two-probe Goda reflection analysis")
    #     plt.grid(True)
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.show()
        
    # ------------------------------------------------------------
    # Enhanced, publication-quality plotting
    # ------------------------------------------------------------
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

        # Avoid the DC point in plots (often zero / noisy)
        flim = f[1:]
        Si_plot = Si[1:]
        Sr_plot = Sr[1:]

        fig, ax = plt.subplots(figsize=(7.0, 3.2), constrained_layout=True)

        ax.plot(
            flim, Si_plot,
            label="Incident spectrum",
            color="#1f77b4",
            linewidth=1.2,
        )
        ax.plot(
            flim, Sr_plot,
            label="Reflected spectrum",
            color="#d62728",
            linestyle="--",
            linewidth=1.2,
        )

        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(r"$S(f)$ [m$^2$/Hz]")
        ax.set_title("Two-probe Goda reflection spectra")
        ax.grid(alpha=0.3)
        ax.legend()

        # Clean look: remove top/right spines
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        # Save in journal-ready formats
        fig.savefig("figures/twoprobe_spectra.png", dpi=300, bbox_inches="tight")
        # fig.savefig("twoprobe_spectra.pdf", bbox_inches="tight")
        plt.show()
        # plt.close(fig)

 
    # Retained energy fraction (after masking) diagnostic
    # Use incident spectrum energy retained vs unmasked incident spectrum energy.
    Si_all = (a_inc ** 2) / (2.0 * df)
    if Si_all.size > 0:
        Si_all = Si_all.copy()
        Si_all[0] = 0.0
    total_energy = np.nansum(Si_all) * df
    retained_energy = np.nansum(Si_int) * df
    retained_energy_fraction = (retained_energy / total_energy) if total_energy > 0 else np.nan

    return {
        "Kr": Kr,
        "Hi": Hi,
        "Hr": Hr,
        "Si": Si,
        "Sr": Sr,
        "f": f,
        "cond": cond_M,
        "bad_cond": bad_cond,
        "retained_energy_fraction": retained_energy_fraction,
    }
