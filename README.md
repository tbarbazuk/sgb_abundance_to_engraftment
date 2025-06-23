# SGB Engraftment Analysis

This repository provides tools for analyzing Species Genome Bin (SGB) engraftment in recipients following Fecal Microbiota Transplantation (FMT), and for associating SGB engraftment frequency with microbial relative abundance measured by MetaPhlAn.

## Features

- Filter metadata by FMT timepoint (e.g., `t7`, `t30`).
- Map recipient samples to donor samples.
- Compute engrafted SGB counts and fractions per donor.
- Export engraftment metrics for downstream analysis.
- Correlate SGB engraftment frequency with MetaPhlAn abundance estimates.
- Command-line argument support for timepoint flexibility.

## Requirements

- Python 3.8+
- Required packages:
  - `pandas`
  - `numpy`

You can install dependencies with:

```bash
pip install -r requirements.txt
