"""sensitivity.py

Utilities to *quantify* spacing sensitivity for two-probe and three-probe
incident–reflected decomposition.

These helpers were added to directly address reviewer requests for a more
systematic assessment of cases where:
  - only one probe pair meets the Goda spacing guideline,
  - two probe pairs fail but one remains marginally acceptable,
  - all probe pairs violate the spacing requirement.

The functions here generate synthetic (known-truth) gauge records using
linear wave theory and then run WaveLabX algorithms to compute errors.

Notes
-----
* These are research/validation utilities, not required for routine use.
* They are intentionally lightweight (NumPy only + WaveLabX itself).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .core import compute_wavelength
from .two_probe import two_probe_goda
from .three_probe import three_probe_array


@dataclass(frozen=True)
class SyntheticTruth:
    """Ground-truth parameters used to generate synthetic gauges."""

    Hi: float
    Hr: float
    Kr: float
    fp: float
    h: float
    gpos: tuple[float, float, float]


def _synthetic_irregular_gauges(
    *,
    fs: float,
    duration: float,
    h: float,
    gpos: tuple[float, float, float],
    Tpeak: float,
    Hi: float,
    Kr: float,
    n_comp: int = 64,
    rel_bandwidth: float = 0.25,
    seed: int = 123,
    noise_std: float = 0.0,
) -> tuple[np.ndarray, SyntheticTruth]:
    """Generate 3-gauge synthetic time series with known incident/reflected energy.

    The signal is constructed as a narrowband random-phase superposition around
    the peak frequency. A reflected component is added with amplitude ratio Kr.
    """

    rng = np.random.default_rng(seed)
    fs = float(fs)
    N = int(round(duration * fs))
    t = np.arange(N) / fs

    fp = 1.0 / float(Tpeak)
    fmin = max(0.01, fp * (1.0 - rel_bandwidth))
    fmax = fp * (1.0 + rel_bandwidth)
    f = np.linspace(fmin, fmax, n_comp)

    # Convert target Hi and Kr into a set of component amplitudes.
    # For a sum of cosines with random phases, variance ~= sum(a_i^2)/2.
    # We pick equal-energy components to match the requested Hi.
    Hi = float(Hi)
    Kr = float(Kr)
    Hr = Kr * Hi

    Ei = (Hi / 4.0) ** 2
    Er = (Hr / 4.0) ** 2
    # Equal-energy components
    ai = np.sqrt(2.0 * Ei / n_comp) * np.ones_like(f)
    ar = np.sqrt(2.0 * Er / n_comp) * np.ones_like(f)

    # Random phases
    ph_i = rng.uniform(0, 2 * np.pi, size=n_comp)
    ph_r = rng.uniform(0, 2 * np.pi, size=n_comp)

    # Wavenumber using linear dispersion (via wavelength)
    k = np.zeros_like(f)
    for idx, fi in enumerate(f):
        L = compute_wavelength(h, 1.0 / fi)
        k[idx] = 2.0 * np.pi / L

    eta = np.zeros((N, 3), dtype=float)
    for gi, x in enumerate(gpos):
        # incident: cos(ωt - kx + φ)
        # reflected: cos(ωt + kx + φ)
        sig = np.zeros(N)
        for ii, fi in enumerate(f):
            w = 2.0 * np.pi * fi
            sig += ai[ii] * np.cos(w * t - k[ii] * x + ph_i[ii])
            sig += ar[ii] * np.cos(w * t + k[ii] * x + ph_r[ii])
        if noise_std > 0:
            sig += rng.normal(0.0, noise_std, size=N)
        eta[:, gi] = sig

    truth = SyntheticTruth(
        Hi=Hi,
        Hr=Hr,
        Kr=(Hr / Hi) if Hi > 0 else np.nan,
        fp=fp,
        h=h,
        gpos=gpos,
    )
    return eta, truth


def spacing_sensitivity(
    *,
    fs: float = 100.0,
    duration: float = 200.0,
    h: float = 0.25,
    Tpeak: float = 1.33,
    Hi: float = 0.035,
    Kr: float = 0.15,
    gpos_sets: Iterable[tuple[float, float, float]] | None = None,
    seed: int = 123,
    noise_std: float = 0.0,
    min_retained_energy: float = 0.8,
) -> dict:
    """Run a spacing sensitivity sweep using synthetic gauges.

    Parameters
    ----------
    gpos_sets:
        Iterable of 3-tuples (x1, x2, x3) in meters.
        If None, a small default set is used that covers:
          - one pair admissible,
          - marginal pairs,
          - all pairs inadmissible.

    Returns
    -------
    dict with:
      - 'truth': ground-truth settings
      - 'results': list of per-geometry results (two-probe per pair + three-probe)
    """

    if gpos_sets is None:
        # These are illustrative defaults; users should customize for their lab.
        # x1 is set to 0.0; x2 and x3 vary.
        gpos_sets = [
            (0.0, 0.15, 0.30),   # typically good (two pairs often admissible)
            (0.0, 0.05, 0.10),   # very small spacing (often violates 0.05L)
            (0.0, 0.35, 0.70),   # large spacing (can violate 0.45L)
            (0.0, 0.10, 0.60),   # one pair admissible, others may not
            (0.0, 0.60, 0.90),   # the paper example
        ]

    results = []

    # Representative wavelength for reporting dx/Lp
    Lp = compute_wavelength(h, Tpeak)

    for idx, gpos in enumerate(gpos_sets):
        eta, truth = _synthetic_irregular_gauges(
            fs=fs,
            duration=duration,
            h=h,
            gpos=gpos,
            Tpeak=Tpeak,
            Hi=Hi,
            Kr=Kr,
            seed=seed + idx,
            noise_std=noise_std,
        )

        # Two-probe: evaluate all 3 pairs
        pairs = [(0, 1), (0, 2), (1, 2)]
        two = []
        for (i, j) in pairs:
            dx = abs(gpos[j] - gpos[i])
            dx_over_Lp = dx / Lp if Lp > 0 else np.nan
            out = two_probe_goda(
                eta[:, [i, j]],
                fs=fs,
                h=h,
                gpos=(gpos[i], gpos[j]),
                plot=False,
                window=None,
            )
            out.update(
                {
                    "pair": (i + 1, j + 1),
                    "dx": dx,
                    "dx_over_Lp": float(dx_over_Lp),
                    "goda_admissible": bool(0.05 <= dx_over_Lp <= 0.45),
                    "err_Hi": float(out["Hi"] - truth.Hi),
                    "err_Hr": float(out["Hr"] - truth.Hr),
                    "err_Kr": float(out["Kr"] - truth.Kr),
                }
            )
            two.append(out)

        # Three-probe
        th = three_probe_array(
            eta,
            fs=fs,
            h=h,
            gpos=gpos,
            plot=False,
            window=None,
            min_retained_energy=min_retained_energy,
        )
        th.update(
            {
                "err_Hi": float(th["Hi"] - truth.Hi),
                "err_Hr": float(th["Hr"] - truth.Hr),
                "err_Kr": float(th["Kr"] - truth.Kr),
            }
        )

        results.append(
            {
                "gpos": gpos,
                "Lp": float(Lp),
                "two_probe": two,
                "three_probe": th,
            }
        )

    return {
        "truth": {
            "Hi": Hi,
            "Kr": Kr,
            "Hr": Kr * Hi,
            "Tpeak": Tpeak,
            "Lp": float(Lp),
            "fs": fs,
            "duration": duration,
            "h": h,
            "noise_std": noise_std,
        },
        "results": results,
    }
