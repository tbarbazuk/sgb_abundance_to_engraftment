import argparse
import pandas as pd
import re
from scipy.stats import spearmanr

def extract_sgb_number(s):
    """Extract digits from an SGB string, return as string."""
    match = re.search(r'(\d+)', str(s))
    return match.group(1) if match else None

def main(engraftment_csv, metaphlan_tsv, donor_prefix, output_csv):
    # Load engraftment CSV
    engraftment_df = pd.read_csv(engraftment_csv)
    engraftment_df['sgb_num'] = engraftment_df['sgb_id'].apply(extract_sgb_number)

    # Load MetaPhlAn TSV - read header rows to get sample names
    print("Reading MetaPhlAn file header...")
    with open(metaphlan_tsv, 'r') as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()
    
    # Parse the second line to get actual sample names
    sample_names = second_line.split('\t')[1:]  # Skip first column (SGB names)
    print(f"Found {len(sample_names)} total samples in MetaPhlAn file")
    
    # Find donor sample indices and names
    donor_indices = []
    donor_names = []
    for i, name in enumerate(sample_names):
        if isinstance(name, str) and name.startswith(donor_prefix):
            donor_indices.append(i + 1)  # +1 because we skip the first column (SGB names)
            donor_names.append(name)
    
    if not donor_names:
        print(f"ERROR: No samples found with prefix '{donor_prefix}'")
        print(f"Available samples: {sample_names[:10]}...")  # Show first 10 samples
        return
    
    print(f"Found {len(donor_names)} donor samples: {donor_names}")

    # Load the full MetaPhlAn file, skipping the first two rows
    print("Loading MetaPhlAn abundance data...")
    metaphlan_df = pd.read_csv(metaphlan_tsv, sep="\t", header=None, skiprows=2, low_memory=False)
    
    # Extract SGB names (first column)
    sgb_names = metaphlan_df.iloc[:, 0].reset_index(drop=True)
    sgb_nums = sgb_names.apply(extract_sgb_number)
    
    # Extract donor columns - need to be careful with indexing
    donor_ra_matrix = metaphlan_df.iloc[:, donor_indices].copy()
    donor_ra_matrix.columns = donor_names
    
    # Convert to numeric, handling any non-numeric values
    donor_ra_matrix = donor_ra_matrix.apply(pd.to_numeric, errors='coerce')
    
    # Add SGB identifier columns
    donor_ra_matrix.insert(0, 'sgb_num', sgb_nums)
    donor_ra_matrix.insert(1, 'sgb_id', sgb_names)
    
    print("\n Verifying donor relative abundances (first 5 SGBs, all donor samples):")
    print(donor_ra_matrix.iloc[:5])

    # Compute mean RA across donor samples
    donor_ra_matrix['relative_abundance'] = donor_ra_matrix[donor_names].mean(axis=1, skipna=True)

    print(f"\n Summary statistics for donor relative abundances:")
    print(f"   • Total SGBs: {len(donor_ra_matrix)}")
    print(f"   • SGBs with mean RA > 0: {(donor_ra_matrix['relative_abundance'] > 0).sum()}")
    print(f"   • Mean RA range: {donor_ra_matrix['relative_abundance'].min():.6f} - {donor_ra_matrix['relative_abundance'].max():.6f}")

    # Debug: Check for inconsistencies
    suspicious = donor_ra_matrix[
        (donor_ra_matrix['relative_abundance'] == 0) &
        (donor_ra_matrix[donor_names].sum(axis=1) > 0)
    ]
    if not suspicious.empty:
        print(f"\n Warning: {len(suspicious)} SGBs have mean RA == 0 but non-zero individual donor RAs")

    # Remove rows with all NaN values and create final donor abundance DataFrame
    donor_abundance_df = donor_ra_matrix[['sgb_id', 'sgb_num', 'relative_abundance']].dropna()
    print(f" Final donor abundance dataset: {len(donor_abundance_df)} SGBs")

    # Identify top 10 abundant SGBs in the donor
    top10_df = donor_abundance_df.sort_values(by="relative_abundance", ascending=False).head(10)
    print(f"\n Top 10 most abundant SGBs in donor:")
    for i, row in top10_df.iterrows():
        print(f"   {row['sgb_id']}: {row['relative_abundance']:.6f}")

    # Merge top 10 with engraftment data
    top10_with_engraftment = pd.merge(top10_df, engraftment_df, on="sgb_num", how="left")
    top10_with_engraftment["engraftment_fraction_global"] = top10_with_engraftment["engraftment_fraction_global"].fillna(0)

    # Save top 10
    top10_output_csv = output_csv.replace(".csv", "_top10_donor_sgbs.csv")
    top10_with_engraftment.drop(columns=["sgb_num"]).to_csv(top10_output_csv, index=False)
    print(f"\n Top 10 abundant donor SGBs saved to {top10_output_csv}")

    # CRITICAL ANALYSIS: SGBs with zero donor abundance but engraftment
    zero_abundance_engrafted = pd.merge(
        engraftment_df[engraftment_df['engraftment_fraction_global'] > 0], 
        donor_abundance_df, 
        on='sgb_num', 
        how='left'
    )
    # Find SGBs that engraft but have no corresponding donor abundance data
    mysterious_sgbs = zero_abundance_engrafted[zero_abundance_engrafted['relative_abundance'].isna()]
    
    if not mysterious_sgbs.empty:
        print(f"\n ALERT: {len(mysterious_sgbs)} SGBs have engraftment but NO donor abundance data:")
        for _, row in mysterious_sgbs.head(10).iterrows():
            # Use sgb_id from engraftment data, not donor data (which is missing)
            sgb_display = row['sgb_id_x'] if 'sgb_id_x' in row else f"SGB_{row['sgb_num']}"
            print(f"   • {sgb_display}: engraftment = {row['engraftment_fraction_global']:.4f}")
        if len(mysterious_sgbs) > 10:
            print(f"   ... and {len(mysterious_sgbs) - 10} more")
        
        # Save mysterious SGBs for investigation
        mysterious_output = output_csv.replace(".csv", "_mysterious_engraftment.csv")
        mysterious_sgbs_clean = mysterious_sgbs.drop(columns=['sgb_num'], errors='ignore')
        mysterious_sgbs_clean.to_csv(mysterious_output, index=False)
        print(f" Mysterious SGBs saved to {mysterious_output}")
    else:
        print("\n All engrafted SGBs have corresponding donor abundance data")
    
    # Also check for SGBs with very low donor abundance but high engraftment
    if not donor_abundance_df.empty:
        very_low_abundance = donor_abundance_df[donor_abundance_df['relative_abundance'] < 1e-5]  # < 0.001%
        if not very_low_abundance.empty:
            low_ab_high_eng = pd.merge(very_low_abundance, engraftment_df, on='sgb_num', how='inner')
            high_engraftment = low_ab_high_eng[low_ab_high_eng['engraftment_fraction_global'] > 0.1]  # > 10% engraftment
            
            if not high_engraftment.empty:
                print(f"\n INTERESTING: {len(high_engraftment)} SGBs with very low donor abundance (<0.001%) but high engraftment (>10%):")
                for _, row in high_engraftment.head(5).iterrows():
                    # Handle potential column name conflicts from merge
                    sgb_display = row.get('sgb_id_x', row.get('sgb_id_y', row.get('sgb_id', f"SGB_{row['sgb_num']}")))
                    print(f"   • {sgb_display}: abundance = {row['relative_abundance']:.2e}, engraftment = {row['engraftment_fraction_global']:.3f}")
                
                # Save for further investigation
                outlier_output = output_csv.replace(".csv", "_low_abundance_high_engraftment.csv")
                high_engraftment_clean = high_engraftment.drop(columns=['sgb_num'], errors='ignore')
                high_engraftment_clean.to_csv(outlier_output, index=False)
                print(f" Low abundance, high engraftment SGBs saved to {outlier_output}")

    # Merge for matched (engrafted only)
    matched_df = pd.merge(engraftment_df, donor_abundance_df, on='sgb_num', how='inner')
    print(f"\n Matched SGBs (present in both datasets): {len(matched_df)}")
    
    if len(matched_df) > 1:
        spearman_res = spearmanr(matched_df['engraftment_fraction_global'], matched_df['relative_abundance'])
        print(f" Spearman correlation: r = {spearman_res.correlation:.3f}, p = {spearman_res.pvalue:.3e}")
    else:
        print(" Too few matched SGBs for correlation analysis")

    matched_df.drop(columns=['sgb_num'], inplace=True)
    matched_df.to_csv(output_csv, index=False)
    print(f"Merged data (engrafted only) saved to {output_csv}")

    # Merge all donor SGBs with engraftment data (fill 0 for missing)
    all_df = pd.merge(donor_abundance_df, engraftment_df, on='sgb_num', how='left')
    all_df['engraftment_fraction_global'] = all_df['engraftment_fraction_global'].fillna(0)
    all_df.drop(columns=['sgb_num'], inplace=True)

    # Save full merged output
    full_output_csv = output_csv.replace(".csv", "_with_nonengrafted.csv")
    all_df.to_csv(full_output_csv, index=False)
    print(f" Full data (including non-engrafted SGBs) saved to {full_output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Associate engraftment counts with MetaPhlAn relative abundances")
    parser.add_argument('--engraftment_csv', required=True, help="CSV file with engraftment frequencies")
    parser.add_argument('--metaphlan_tsv', required=True, help="MetaPhlAn TSV file")
    parser.add_argument('--output_csv', required=True, help="Output CSV filename for merged data (engrafted only)")
    parser.add_argument('--donor_prefix', required=True, help="Prefix to match donor samples (e.g., 'FMT_DTT001')")

    args = parser.parse_args()
    main(args.engraftment_csv, args.metaphlan_tsv, args.donor_prefix, args.output_csv)
