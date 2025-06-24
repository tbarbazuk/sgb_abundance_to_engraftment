# 🧬 SGB Engraftment Analysis

This repository provides a suite of Python tools to analyze **Species Genome Bin (SGB) engraftment** in recipients following **Fecal Microbiota Transplantation (FMT)**. It includes utilities for parsing sample metadata, identifying engrafted strains, and correlating engraftment patterns with taxonomic profiles from **MetaPhlAn**.

---

## Features

 **Metadata filtering by FMT timepoint** (`t7`, `t30`, etc.)
 **Donor–recipient sample mapping** using sample sheet metadata
 **Engraftment analysis**:

  * Count and fraction of engrafted SGBs per donor
  * Output of matched strain tables
  *  **Correlation of SGB engraftment frequency** with **MetaPhlAn** relative abundance profiles
  *  **Command-line interface** with timepoint argument support for flexible analysis
  *  Structured for modularity and reuse across pipelines

---

## 📂 Script Overview

| Script Name                       | Description                                                                                           |
| --------------------------------- | ----------------------------------------------------------------------------------------------------- |
| `sgb_engraftment.py`              | Core logic for calculating engrafted SGB counts and fractions, and matching strains to donor profiles |
| `engraftment_vs_abundance.py`     | Correlates SGB engraftment frequency with MetaPhlAn taxonomic abundance                               |
| `filter_metadata_by_timepoint.py` | Utility to filter sample metadata for a specified timepoint (`t7`, `t30`, etc.)                       |
| `common.py` *(optional)*          | (Recommended) Shared utility functions like strain parsing, sample ID logic, etc. (if modularized)    |

---

## Inputs

* `metadata.csv`: Sample metadata with columns like `sample_id`, `subject_id`, `donor_sample_name`, and timepoint.
* `engraftment_table.csv`: Matrix or long-format table containing `sample_id`, `eng_strains`, and possibly donor info.
* `metaphlan.tsv` (optional): MetaPhlAn taxonomic relative abundance table.

---

## Outputs

* `engraftment_metrics.csv`: Number and fraction of engrafted SGBs per donor sample.
* `matched_strains.csv`: Table listing engrafted strains present in both donor and recipient.
* `engraftment_vs_abundance.tsv`: (Optional) Correlation statistics for SGB engraftment frequency vs. MetaPhlAn abundance.

---

## 🛠️ Requirements

* Python **3.8+**
* Python packages:

  * `pandas`
  * `numpy`

---

## ⚙️ Example Usage

```bash
# Run engraftment analysis for timepoint t7
python sgb_engraftment.py --timepoint t7 --metadata metadata.csv --engraftment engraftment_table.csv

# Correlate engraftment frequency with MetaPhlAn abundance
python engraftment_vs_abundance.py --engraftment engraftment_metrics.csv --metaphlan metaphlan.tsv
```

