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

    # Extract numeric part from engraftment sgb_id for matching
    engraftment_df['sgb_num'] = engraftment_df['sgb_id'].apply(extract_sgb_number)

    metaphlan_df = pd.read_csv(metaphlan_tsv, sep="\t", header=None, skiprows=1, low_memory=False)
    sample_names = metaphlan_df.iloc[0, 1:]
    print("MetaPhlAn sample names row:")
    print(sample_names.values)

    if donor_sample not in sample_names.values:
        raise ValueError(f"Sample '{donor_sample}' not found in MetaPhlAn row 2.")

    # Find the donor sample column index (relative to metaphlan_df)
    donor_col_idx = sample_names[sample_names == donor_sample].index[0]

    # Extract SGB names from first column starting from row 2 (index 2 after skiprows)
    sgb_names = metaphlan_df.iloc[2:, 0].reset_index(drop=True)
    # Extract numeric portion of sgb names for matching
    sgb_nums = sgb_names.apply(extract_sgb_number)

    # Extract relative abundances for donor sample column (starting from row 2)
    rel_abunds = metaphlan_df.iloc[2:, donor_col_idx].reset_index(drop=True)
    rel_abunds = pd.to_numeric(rel_abunds, errors='coerce')

    # Create DataFrame for donor abundance with numeric sgb ids
    donor_abundance_df = pd.DataFrame({
        'sgb_id': sgb_names,
        'sgb_num': sgb_nums,
        'relative_abundance': rel_abunds
    }).dropna(subset=['relative_abundance', 'sgb_num'])

    print("MetaPhlAn SGB numeric IDs preview:")
    print(donor_abundance_df[['sgb_id', 'sgb_num']].head(10))

    # Merge engraftment data with donor abundance on numeric sgb id
    merged_df = pd.merge(engraftment_df, donor_abundance_df, on='sgb_num', how='inner')

    print(f"Number of matched SGBs: {len(merged_df)}")

    # Calculate Spearman correlation between engraftment fraction and relative abundance
    spearman_res = spearmanr(merged_df['engraftment_fraction_global'], merged_df['relative_abundance'])
    print(f"Spearman correlation: r = {spearman_res.correlation:.3f}, p = {spearman_res.pvalue:.3e}")

    # Drop helper columns and save merged data
    merged_df.drop(columns=['sgb_num'], inplace=True)
    merged_df.to_csv(output_csv, index=False)
    print(f"Merged data saved to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Associate engraftment counts with MetaPhlAn relative abundances")
    parser.add_argument('--engraftment_csv', required=True, help="CSV file with engraftment frequencies")
    parser.add_argument('--metaphlan_tsv', required=True, help="MetaPhlAn TSV file")
    parser.add_argument('--donor_sample', required=True, help="Donor sample name as in MetaPhlAn row 2")
    parser.add_argument('--output_csv', required=True, help="Output CSV filename for merged data")

    args = parser.parse_args()
    main(args.engraftment_csv, args.metaphlan_tsv, args.donor_sample, args.output_csv)

