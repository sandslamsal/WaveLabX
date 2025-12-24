"""
analysis.py

High-level helpers that combine WaveLabX modules into a practical workflow.
"""

from __future__ import annotations

import numpy as np

from .core import compute_wavelength
from .stats import zero_crossing
from .two_probe import two_probe_goda
from .three_probe import three_probe_array


def reflection_analysis(
    eta: np.ndarray,
    fs: float,
    h: float,
    gpos: tuple[float, ...],
    prefer_three_probe: bool = True,
    min_retained_energy: float = 0.8,
    window: str | None = None,
    min_two_probe_retained_energy: float = 0.8,
) -> dict:
    """
    High-level reflection analysis with diagnostics and automatic method selection.

    What this does (in a reproducible, code-level way):
      1) Runs zero-crossing on probe 1 to compute a representative Tmean and Lp.
      2) Evaluates all 2-probe pairs (when available):
           - checks Δx/Lp guideline,
           - computes conditioning diagnostics (cond array),
           - reports retained spectral energy fraction after masking.
      3) If 3 probes are available, runs the 3-probe algorithm and reports:
           - retained_energy_fraction (how much composite energy remains usable),
           - conditioning flags (bad_cond_pair).
      4) Selects a recommended method:
           - prefer 3-probe if it retains enough energy (>= min_retained_energy)
             and prefer_three_probe=True;
           - otherwise select the best admissible 2-probe pair (highest retained
             energy, then lowest worst-condition number).

    Returns a dictionary with keys:
      - 'Hs', 'Tmean', 'Lp' from zero-crossing (probe 1)
      - 'two_probe' (optional), 'three_probe' (optional)
      - 'method_used'
    """
    eta = np.asarray(eta, dtype=float)
    if eta.ndim != 2:
        raise ValueError("eta must be a 2D array with shape (N, n_gauges).")

    n_g = eta.shape[1]
    if len(gpos) != n_g:
        raise ValueError("gpos length must match number of columns in eta.")

    # --- Zero-crossing on probe 1 for a representative period
    zc = zero_crossing(eta[:, 0], fs)
    Hs = float(zc["Hs"])
    Tmean = float(zc["Tmean"])
    Lp = float(compute_wavelength(h, Tmean))

    out = {"Hs": Hs, "Tmean": Tmean, "Lp": Lp, "zero_crossing": zc}

    # --- Candidate two-probe pairs (evaluate all)
    two_probe_results: list[dict] = []
    if n_g >= 2 and Lp > 0:
        for i in range(n_g):
            for j in range(i + 1, n_g):
                dx = abs(gpos[j] - gpos[i])
                r = dx / Lp
                tp = two_probe_goda(
                    eta[:, [i, j]],
                    fs=fs,
                    h=h,
                    gpos=(gpos[i], gpos[j]),
                    plot=False,
                    window=window,
                )
                tp["pair"] = (i + 1, j + 1)  # 1-based
                tp["dx_over_Lp"] = float(r)
                tp["goda_admissible"] = bool(0.05 <= r <= 0.45)
                # a simple scalar conditioning summary
                cond = tp.get("cond", None)
                tp["cond_max"] = float(np.nanmax(cond)) if cond is not None else float("nan")
                two_probe_results.append(tp)

        out["two_probe_all"] = two_probe_results

        # pick the best admissible 2-probe pair
        admissible = [tp for tp in two_probe_results if tp.get("goda_admissible", False)]
        if admissible:
            admissible_sorted = sorted(
                admissible,
                key=lambda d: (
                    -(d.get("retained_energy_fraction", float("nan")) if d.get("retained_energy_fraction", float("nan")) == d.get("retained_energy_fraction", float("nan")) else -1.0),
                    d.get("cond_max", float("inf")),
                    abs(d.get("dx_over_Lp", 9.9) - 0.15),
                ),
            )
            out["two_probe_best"] = admissible_sorted[0]

    # --- Three-probe method (if available)
    if n_g == 3:
        th = three_probe_array(
            eta[:, :3],
            fs=fs,
            h=h,
            gpos=(gpos[0], gpos[1], gpos[2]),
            plot=False,
            window=window,
            min_retained_energy=min_retained_energy,
        )
        out["three_probe"] = th

        # Automatic method selection
        use_three = False
        if prefer_three_probe:
            frac = th.get("retained_energy_fraction", float("nan"))
            if (frac == frac) and (frac >= min_retained_energy):
                use_three = True

        if use_three:
            out["method_used"] = "three_probe"
        else:
            # fall back to best admissible 2-probe (if it also retains enough energy)
            best2 = out.get("two_probe_best", None)
            if best2 is not None:
                frac2 = best2.get("retained_energy_fraction", float("nan"))
                if (frac2 == frac2) and (frac2 >= min_two_probe_retained_energy):
                    out["method_used"] = "two_probe"
                else:
                    out["method_used"] = "three_probe"  # still return 3-probe result but flagged by diagnostic
                    out["warnings"] = out.get("warnings", []) + [
                        f"Best 2-probe pair retained only {frac2:.2%} energy (< {min_two_probe_retained_energy:.0%}).",
                    ]
            else:
                out["method_used"] = "three_probe"

    else:
        out["method_used"] = "two_probe" if out.get("two_probe_best") is not None else "none"

    return out
