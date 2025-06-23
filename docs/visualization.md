# Visualizing Engraftment vs. Relative Abundance

This script visualizes the relationship between **engraftment frequency** and **relative abundance** using scatter plots, histograms, and log-transformed plots. It also calculates a Spearman correlation coefficient to assess monotonic association between the two variables.

---

## 🔧 Input Requirements

1. **Merged CSV File**  
   - A CSV containing both engraftment frequency and donor relative abundance for each SGB.
   - Typically generated using `relative_abundance_to_engraftment.py`.

2. **Output Prefix**  
   - A filename prefix used to name all output plot files.

---

## ▶️ How to Run

Run the script from the command line as follows:

```bash
python plot_engraftment_relationship.py \
  --merged_csv donor01_merged.csv \
  --output_prefix donor01

