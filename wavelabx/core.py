
"""
core.py

Shared constants and core wave mechanics utilities for WaveLabX.
"""

import numpy as np

GRAVITY = 9.81


def compute_wavelength(h: float, T: float, tol: float = 1e-5, max_iter: int = 50) -> float:
    """Compute linear wave length L for given water depth h and period T
    by solving the dispersion relation iteratively:

        ω² = g k tanh(k h),   where ω = 2π / T,  k = 2π / L

    Parameters
    ----------
    h : float
        Water depth [m].
    T : float
        Wave period [s].
    tol : float, optional
        Convergence tolerance on L [m].
    max_iter : int, optional
        Maximum number of iterations.

    Returns
    -------
    L : float
        Wave length [m].
    """
    g = GRAVITY
    w = 2.0 * np.pi / T

    # Deep-water estimate as initial guess
    L0 = g * T**2 / (2.0 * np.pi)
    L = L0

    for _ in range(max_iter):
        k = 2.0 * np.pi / L
        L_next = (g * T**2 / (2.0 * np.pi)) * np.tanh(k * h)
        if abs(L_next - L) < tol:
            L = L_next
            break
        L = L_next

    return float(L)
