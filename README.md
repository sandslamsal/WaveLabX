# WaveLabX
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18049227.svg)](https://doi.org/10.5281/zenodo.18049227)

# WaveLabX

WaveLabX is an open-source Python toolkit for laboratory wave-probe analysis, providing reproducible wave statistics and incident–reflected decomposition.

- Zero-crossing wave statistics from single-probe records
- Two-probe Goda–Suzuki frequency-domain decomposition
- Three-probe redundant-array decomposition with validity filtering


![WaveLabX workflow and architecture](figures/wavelabx_architecture.jpg)


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


## Browser tool

**Live demo: [wave-lab-x.vercel.app](https://wave-lab-x.vercel.app)**

`web/` contains a self-contained browser application for three-probe
reflection analysis — drop one or more 6-channel wave-gauge CSV files and
get a complete table of incident/reflected wave heights and reflection
coefficients for two probe arrays. It runs entirely client-side (no
install, no server) and is deployable as a static site on Vercel.

To run it locally, open `web/index.html`, or serve the folder:

```bash
cd web && python3 -m http.server 8000
```

To deploy on Vercel, import this repository and set the project
**Root Directory** to `web`.


## Files included

- `wavelabx/` — package source (API documented in docstrings)
- `web/` — browser-based three-probe reflection analysis tool
- `run_wavelabx_example.ipynb` — example notebook
- `wavedata.csv` — example dataset
- `figures/` — figures for manuscript and examples
- `paper.md`, `paper.bib` — JOSS manuscript source


## License & Citation

WaveLabX is released under the MIT License (see `LICENSE`).

If you use WaveLabX in published research, please cite the associated JOSS manuscript and the archived software release on Zenodo. Citation metadata is provided in `CITATION.cff`.

Zenodo archive DOI: https://doi.org/10.5281/zenodo.18049227

