# Engraftment Counter

This script analyzes engraftment data for fecal microbiota transplantation (FMT) samples at specified timepoints.

## Overview

Given a metadata TSV and an engraftment CSV file, the script:
- Filters samples by a user-specified timepoint (e.g., `t7`)
- Maps samples to donors
- Loads engraftment strain data
- Matches strains between donor sample and recipient
- Counts matched strains (engrafted SGBs)
- Computes per-donor and global engraftment fractions
- Outputs summary statistics and a CSV of global engraftment frequencies

## Usage

```bash
python engraftment_counter.py <timepoint> <metadata_file> <engraftment_file>

