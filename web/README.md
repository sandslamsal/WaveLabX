# Wave Reflection Analysis

**Live demo: [wave-lab-x.vercel.app](https://wave-lab-x.vercel.app)**

A browser-based tool for **three-probe wave reflection analysis** —
separation of incident and reflected wave heights from co-linear
capacitance wave-gauge records.

Drop one or more 6-channel time-series CSV files and the tool returns a
complete table of incident/reflected wave heights and reflection
coefficients for two probe arrays.

## What it computes

For each file (6 columns of wave-gauge time series):

| Quantity | From |
|---|---|
| `Hi1`, `Hr1`, `Kr1` | channels 1–3 — seaward probe array |
| `Hi2`, `Hr2`, `Kr2` | channels 4–6 — shoreward probe array |
| `Kt` | `Hi2 / Hi1` — transmission coefficient |
| `f`, `T` | wave frequency / period, detected from the signal |
| `L` | local wavelength from linear dispersion |

## Method

Three-wave-probe separation following **Goda & Suzuki (1976)** and
**Kobayashi, Cox & Wurjanto (1990)**, as applied in:

> Lamsal, S., Haus, B. K. & Rhode-Barbarigos, L. (2026).
> *An experimental study on wave transmission over submerged SEAHIVE® breakwaters.*
> Coastal Engineering Journal.
> [doi:10.1080/21664250.2026.2661171](https://doi.org/10.1080/21664250.2026.2661171)

The complex wave amplitude at the wave frequency is obtained from a
Hann-windowed discrete Fourier transform of each gauge record (DC offset
removed). The wave number follows from the linear dispersion relation
`ω² = g k tanh(k d)`, and incident/reflected amplitudes are recovered by
complex least squares over the three gauges. The method is valid for
gauge spacing `0.05 ≤ Δl/L ≤ 0.45`; outside this range the tool falls
back to the best-conditioned gauge pair.

## Features

- **Drag-and-drop** batch processing of multiple files
- **Automatic wave-frequency detection** from the signal spectrum
  (does not rely on file names)
- **Global water depth** with one-click apply-to-all, plus per-file override
- **Configurable gauge positions** for both probe arrays
- **CSV export** of the full results table
- Runs **entirely client-side** — no data leaves the browser

## Input format

Each CSV file: 6 numeric columns (one per wave gauge), optional header
row. Water depth is read from the file name if it contains
`Depth=<value>`, otherwise the global depth setting is used.

## Running locally

No build step. Open `index.html` directly, or serve the folder:

```bash
python3 -m http.server 8000
# then open http://localhost:8000
```

## Deployment

Static site — deploys on Vercel with no configuration.
