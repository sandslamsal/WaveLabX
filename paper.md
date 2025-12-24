---
title: 'WaveLabX: An open-source Python toolkit for wave-probe statistics and incident–reflected decomposition'
tags:
  - Python
  - coastal engineering
  - wave flume
  - reflection analysis
  - research software
authors:
  - name: Sandesh Lamsal
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Landolf Rhode-Barbarigos
    affiliation: 1
affiliations:
  - name: Department of Civil and Architectural Engineering, University of Miami, USA
    index: 1
date: 2025-XX-XX
bibliography: paper.bib
---

## Summary

WaveLabX is an open-source Python toolkit for processing laboratory wave-probe
time series in coastal and hydraulic engineering. The software integrates
classical zero-crossing wave statistics with two-probe (Goda–Suzuki) and
three-probe array methods for incident–reflected wave decomposition into a
single, reproducible workflow. WaveLabX provides automated probe-spacing,
numerical-conditioning, and retained-energy diagnostics to guide method
selection and improve reliability of reflection analysis in wave flumes and
basins.

## Statement of need

Incident–reflected wave decomposition is a standard task in laboratory wave
experiments, yet practical implementations of classical methods are often
fragmented across laboratories as undocumented scripts. This fragmentation
limits reproducibility, obscures numerical limitations related to probe
spacing, and complicates cross-laboratory comparison of results. WaveLabX
addresses this gap by providing a transparent, validated, and reusable
implementation of established reflection-analysis techniques in a modern
scientific Python environment.

## Software description

WaveLabX is implemented as a modular Python package with separate components
for dispersion calculations, wave statistics, two-probe reflection analysis,
and three-probe array decomposition. The software supports configurable
preprocessing, optional windowing, and frequency-domain diagnostics, and
exposes a clean public API for integration into laboratory workflows and
Jupyter notebooks.

## Illustrative example

WaveLabX includes example notebooks demonstrating analysis of three-gauge wave
records, including computation of bulk wave statistics, probe-spacing
diagnostics, and reconstruction of incident and reflected wave spectra. These
examples reproduce published benchmark cases and laboratory measurements with
sub-percent error.

## Quality control

The software includes automated tests for core numerical routines, validation
against analytical and published reference cases, and continuous integration
to ensure reproducibility across environments. Diagnostics such as condition
numbers and retained spectral energy are reported to warn users when classical
reflection methods become unreliable.

## Availability

WaveLabX is released under the MIT License and is available at:
https://github.com/sandslamsal/WaveLabX

The version associated with this paper has been archived on Zenodo and assigned
a DOI.

## References
