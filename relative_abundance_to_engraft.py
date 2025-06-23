import argparse
import pandas as pd
import re
from scipy.stats import spearmanr

def extract_sgb_number(s):
    """Extract digits from an SGB string, return as string."""
    match = re.search(r'(\d+)', str(s))
    return match.group(1) if match else None

def main(engraftment_csv, metaphlan_tsv, donor_sample, output_csv):
    # Load engraftment CSV
    engraftment_df = pd.read_csv(engraftment_csv)
    print("Original engraftment SGB IDs preview:")
    print(engraftment_df['sgb_id'].head(10))

    engraftment_df['sgb_num'] = engraftment_df['sgb_id'].apply(extract_sgb_number)

    metaphlan_df = pd.read_csv(metaphlan_tsv, sep="\t", header=None, skiprows=1, low_memory=False)
    sample_names = metaphlan_df.iloc[0, 1:]
    print("MetaPhlAn sample names row:")
    print(sample_names.values)

    if donor_sample not in sample_names.values:
        raise ValueError(f"Sample '{donor_sample}' not found in MetaPhlAn row 2.")

    donor_col_idx = sample_names[sample_names == donor_sample].index[0]
    sgb_names = metaphlan_df.iloc[2:, 0].reset_index(drop=True)
    sgb_nums = sgb_names.apply(extract_sgb_number)
    rel_abunds = metaphlan_df.iloc[2:, donor_col_idx].reset_index(drop=True)
    rel_abunds = pd.to_numeric(rel_abunds, errors='coerce')

    donor_abundance_df = pd.DataFrame({
        'sgb_id': sgb_names,
        'sgb_num': sgb_nums,
        'relative_abundance': rel_abunds
    }).dropna(subset=['relative_abundance', 'sgb_num'])

    print("MetaPhlAn SGB numeric IDs preview:")
    print(donor_abundance_df[['sgb_id', 'sgb_num']].head(10))

    # 1. Engrafted only (inner join)
    matched_df = pd.merge(engraftment_df, donor_abundance_df, on='sgb_num', how='inner')
    print(f"Number of matched SGBs: {len(matched_df)}")

    spearman_res = spearmanr(matched_df['engraftment_fraction_global'], matched_df['relative_abundance'])
    print(f"Spearman correlation: r = {spearman_res.correlation:.3f}, p = {spearman_res.pvalue:.3e}")

    matched_df.drop(columns=['sgb_num'], inplace=True)
    matched_df.to_csv(output_csv, index=False)
    print(f"Merged data (engrafted only) saved to {output_csv}")

    # 2. All donor SGBs with engraftment if present (left join)
    all_df = pd.merge(donor_abundance_df, engraftment_df, on='sgb_num', how='left')
    all_df['engraftment_fraction_global'] = all_df['engraftment_fraction_global'].fillna(0)
    all_df.drop(columns=['sgb_num'], inplace=True)

    # Save to additional CSV
    full_output_csv = output_csv.replace(".csv", "_with_nonengrafted.csv")
    all_df.to_csv(full_output_csv, index=False)
    print(f"Full data (including non-engrafted SGBs) saved to {full_output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Associate engraftment counts with MetaPhlAn relative abundances")
    parser.add_argument('--engraftment_csv', required=True, help="CSV file with engraftment frequencies")
    parser.add_argument('--metaphlan_tsv', required=True, help="MetaPhlAn TSV file")
    parser.add_argument('--donor_sample', required=True, help="Donor sample name as in MetaPhlAn row 2")
    parser.add_argument('--output_csv', required=True, help="Output CSV filename for merged data (engrafted only)")

    args = parser.parse_args()
    main(args.engraftment_csv, args.metaphlan_tsv, args.donor_sample, args.output_csv)
