
"""
WaveLabX: Wave probe analysis toolkit

Public API:
- compute_wavelength
- zero_crossing
- two_probe_goda
- three_probe_array
"""

from .core import compute_wavelength, GRAVITY
from .stats import zero_crossing
from .two_probe import two_probe_goda
from .three_probe import three_probe_array
from .analysis import reflection_analysis
from .sensitivity import spacing_sensitivity

__all__ = [
    "GRAVITY",
    "compute_wavelength",
    "zero_crossing",
    "two_probe_goda",
    "three_probe_array",
    "reflection_analysis",
    "spacing_sensitivity",
]
