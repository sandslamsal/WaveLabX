# WaveLabX
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18049227.svg)](https://doi.org/10.5281/zenodo.18049227)

# WaveLabX

WaveLabX is an open-source Python toolkit for laboratory wave-probe analysis, providing reproducible wave statistics and incident–reflected decomposition.

- Zero-crossing wave statistics from single-probe records
- Two-probe Goda–Suzuki frequency-domain decomposition
- Three-probe redundant-array decomposition with validity filtering


![WaveLabX workflow and architecture](figures/wavelabx_architecture.pdf)


## Installation

Recommended: use a Python virtual environment.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```


## Minimal usage

```python
import wavelabx
# See run_wavelabx_example.ipynb for a full demonstration
```


## Files included

- `wavelabx/` — package source (API documented in docstrings)
- `run_wavelabx_example.ipynb` — example notebook
- `wavedata.csv` — example dataset
- `figures/` — figures for manuscript and examples
- `paper.md`, `paper.bib` — JOSS manuscript source


## License & Citation

WaveLabX is released under the MIT License (see `LICENSE`).

If you use WaveLabX in published research, please cite the JOSS manuscript (paper.md) and see CITATION.cff for citation details. A Zenodo DOI will be added after the first release.
