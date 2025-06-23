#  Merging Engraftment and MetaPhlAn Abundance Data

This script merges engraftment frequency data (generated using `engraftment_counter.py`) with relative abundance data from MetaPhlAn output. The goal is to associate each SGB's engraftment rate with its abundance in the donor sample, allowing downstream statistical analysis and visualization.

---

##  Input Requirements

1. **Engraftment CSV**  
   - A CSV file containing engraftment frequencies for each SGB.
   - Typically generated using `engraftment_counter.py` for a specific timepoint.

2. **MetaPhlAn TSV**  
   - A MetaPhlAn output `.tsv` file with relative abundance data for each SGB, associated with FMT samples and their respective timepoints.
   - Modifications may be necessary depending on the formatting of your MetaPhlAn output file.

3. **Donor Sample Name**  
   - The sample ID from MetaPhlAn corresponding to the donor.

---

## ▶️ How to Run

From the command line, run:

```bash
python merge_engraftment_and_abundance.py \
  --engraftment_csv engraftment_frequencies_t7.csv \
  --metaphlan_tsv metaphlan_profiles_sgb.tsv \
  --donor_sample FMT_Donor1_merged \
  --output_csv donor01_merged.csv

