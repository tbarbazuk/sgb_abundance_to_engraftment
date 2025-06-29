import argparse
import pandas as pd
import re
from scipy.stats import spearmanr

def extract_sgb_number(s):
    """Extract digits from an SGB string, return as string."""
    match = re.search(r'(\d+)', str(s))
    return match.group(1) if match else None

def load_sample_list(sample_list_file):
    """Load sample names from a text file (one per line)."""
    with open(sample_list_file, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]
    return samples

def main(engraftment_csv, metaphlan_tsv, sample_list_file, output_csv):
    # Load target sample names
    target_samples = load_sample_list(sample_list_file)
    print(f"Target samples from {sample_list_file}: {len(target_samples)} samples")
    print(f"First few samples: {target_samples[:5]}")

    # Load engraftment CSV
    engraftment_df = pd.read_csv(engraftment_csv)
    engraftment_df['sgb_num'] = engraftment_df['sgb_id'].apply(extract_sgb_number)
    
    # Print columns to debug
    print(f"Engraftment CSV columns: {list(engraftment_df.columns)}")

    # Load MetaPhlAn TSV - read header rows to get sample names
    print("Reading MetaPhlAn file header...")
    with open(metaphlan_tsv, 'r') as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()
    
    # Parse the second line to get actual sample names
    sample_names = second_line.split('\t')[1:]  # Skip first column (SGB names)
    print(f"Found {len(sample_names)} total samples in MetaPhlAn file")
    
    # Find target sample indices and names
    target_indices = []
    found_samples = []
    missing_samples = []
    
    for target_sample in target_samples:
        try:
            idx = sample_names.index(target_sample)
            target_indices.append(idx + 1)  # +1 because we skip the first column (SGB names)
            found_samples.append(target_sample)
        except ValueError:
            missing_samples.append(target_sample)
    
    print(f"Found {len(found_samples)}/{len(target_samples)} target samples in MetaPhlAn file")
    
    if missing_samples:
        print(f"WARNING: {len(missing_samples)} samples not found in MetaPhlAn file:")
        for sample in missing_samples[:10]:  # Show first 10 missing
            print(f"   • {sample}")
        if len(missing_samples) > 10:
            print(f"   ... and {len(missing_samples) - 10} more")
    
    if not found_samples:
        print("ERROR: No target samples found in MetaPhlAn file")
        print(f"Available samples: {sample_names[:10]}...")  # Show first 10 samples
        return
    
    print(f"Using {len(found_samples)} samples for analysis: {found_samples}")

    # Load the full MetaPhlAn file, skipping the first two rows
    print("Loading MetaPhlAn abundance data...")
    metaphlan_df = pd.read_csv(metaphlan_tsv, sep="\t", header=None, skiprows=2, low_memory=False)
    
    # Extract SGB names (first column)
    sgb_names = metaphlan_df.iloc[:, 0].reset_index(drop=True)
    sgb_nums = sgb_names.apply(extract_sgb_number)
    
    # Extract target sample columns
    target_ra_matrix = metaphlan_df.iloc[:, target_indices].copy()
    target_ra_matrix.columns = found_samples
    
    # Convert to numeric, handling any non-numeric values
    target_ra_matrix = target_ra_matrix.apply(pd.to_numeric, errors='coerce')
    
    # Add SGB identifier columns
    target_ra_matrix.insert(0, 'sgb_num', sgb_nums)
    target_ra_matrix.insert(1, 'sgb_id', sgb_names)
    
    print("\nVerifying target sample relative abundances (first 5 SGBs, all target samples):")
    print(target_ra_matrix.iloc[:5])

    # Compute mean RA across target samples
    target_ra_matrix['relative_abundance'] = target_ra_matrix[found_samples].mean(axis=1, skipna=True)

    print(f"\nSummary statistics for target sample relative abundances:")
    print(f"   • Total SGBs: {len(target_ra_matrix)}")
    print(f"   • SGBs with mean RA > 0: {(target_ra_matrix['relative_abundance'] > 0).sum()}")
    print(f"   • Mean RA range: {target_ra_matrix['relative_abundance'].min():.6f} - {target_ra_matrix['relative_abundance'].max():.6f}")

    # Debug: Check for inconsistencies
    suspicious = target_ra_matrix[
        (target_ra_matrix['relative_abundance'] == 0) &
        (target_ra_matrix[found_samples].sum(axis=1) > 0)
    ]
    if not suspicious.empty:
        print(f"\nWarning: {len(suspicious)} SGBs have mean RA == 0 but non-zero individual sample RAs")

    # Remove rows with all NaN values and create final abundance DataFrame
    sample_abundance_df = target_ra_matrix[['sgb_id', 'sgb_num', 'relative_abundance']].dropna()
    print(f"Final sample abundance dataset: {len(sample_abundance_df)} SGBs")

    # Identify top 10 abundant SGBs in the target samples
    top10_df = sample_abundance_df.sort_values(by="relative_abundance", ascending=False).head(10)
    print(f"\nTop 10 most abundant SGBs in target samples:")
    for i, row in top10_df.iterrows():
        print(f"   {row['sgb_id']}: {row['relative_abundance']:.6f}")

    # Merge top 10 with engraftment data
    top10_with_engraftment = pd.merge(top10_df, engraftment_df, on="sgb_num", how="left")
    # FIXED: Use correct column name
    top10_with_engraftment["engraftment_frequency"] = top10_with_engraftment["engraftment_frequency"].fillna(0)

    # Save top 10
    top10_output_csv = output_csv.replace(".csv", "_top10_target_sgbs.csv")
    top10_with_engraftment.drop(columns=["sgb_num"]).to_csv(top10_output_csv, index=False)
    print(f"\nTop 10 abundant target sample SGBs saved to {top10_output_csv}")

    # CRITICAL ANALYSIS: SGBs with zero sample abundance but engraftment
    zero_abundance_engrafted = pd.merge(
        engraftment_df[engraftment_df['engraftment_frequency'] > 0], 
        sample_abundance_df, 
        on='sgb_num', 
        how='left'
    )
    # Find SGBs that engraft but have no corresponding sample abundance data
    mysterious_sgbs = zero_abundance_engrafted[zero_abundance_engrafted['relative_abundance'].isna()]
    
    if not mysterious_sgbs.empty:
        print(f"\nALERT: {len(mysterious_sgbs)} SGBs have engraftment but NO sample abundance data:")
        for _, row in mysterious_sgbs.head(10).iterrows():
            # Use sgb_id from engraftment data, not sample data (which is missing)
            sgb_display = row['sgb_id_x'] if 'sgb_id_x' in row else f"SGB_{row['sgb_num']}"
            print(f"   • {sgb_display}: engraftment = {row['engraftment_frequency']:.4f}")
        if len(mysterious_sgbs) > 10:
            print(f"   ... and {len(mysterious_sgbs) - 10} more")
        
        # Save mysterious SGBs for investigation
        mysterious_output = output_csv.replace(".csv", "_mysterious_engraftment.csv")
        mysterious_sgbs_clean = mysterious_sgbs.drop(columns=['sgb_num'], errors='ignore')
        mysterious_sgbs_clean.to_csv(mysterious_output, index=False)
        print(f"Mysterious SGBs saved to {mysterious_output}")
    else:
        print("\nAll engrafted SGBs have corresponding sample abundance data")
    
    # Also check for SGBs with very low sample abundance but high engraftment
    if not sample_abundance_df.empty:
        very_low_abundance = sample_abundance_df[sample_abundance_df['relative_abundance'] < 1e-5]  # < 0.001%
        if not very_low_abundance.empty:
            low_ab_high_eng = pd.merge(very_low_abundance, engraftment_df, on='sgb_num', how='inner')
            high_engraftment = low_ab_high_eng[low_ab_high_eng['engraftment_frequency'] > 0.1]  # > 10% engraftment
            
            if not high_engraftment.empty:
                print(f"\nINTERESTING: {len(high_engraftment)} SGBs with very low sample abundance (<0.001%) but high engraftment (>10%):")
                for _, row in high_engraftment.head(5).iterrows():
                    # Handle potential column name conflicts from merge
                    sgb_display = row.get('sgb_id_x', row.get('sgb_id_y', row.get('sgb_id', f"SGB_{row['sgb_num']}")))
                    print(f"   • {sgb_display}: abundance = {row['relative_abundance']:.2e}, engraftment = {row['engraftment_frequency']:.3f}")
                
                # Save for further investigation
                outlier_output = output_csv.replace(".csv", "_low_abundance_high_engraftment.csv")
                high_engraftment_clean = high_engraftment.drop(columns=['sgb_num'], errors='ignore')
                high_engraftment_clean.to_csv(outlier_output, index=False)
                print(f"Low abundance, high engraftment SGBs saved to {outlier_output}")

    # Merge for matched (engrafted only)
    matched_df = pd.merge(engraftment_df, sample_abundance_df, on='sgb_num', how='inner')
    print(f"\nMatched SGBs (present in both datasets): {len(matched_df)}")
    
    if len(matched_df) > 1:
        spearman_res = spearmanr(matched_df['engraftment_frequency'], matched_df['relative_abundance'])
        print(f"Spearman correlation: r = {spearman_res.correlation:.3f}, p = {spearman_res.pvalue:.3e}")
    else:
        print("Too few matched SGBs for correlation analysis")

    matched_df.drop(columns=['sgb_num'], inplace=True)
    matched_df.to_csv(output_csv, index=False)
    print(f"Merged data (engrafted only) saved to {output_csv}")

    # Merge all sample SGBs with engraftment data (fill 0 for missing)
    all_df = pd.merge(sample_abundance_df, engraftment_df, on='sgb_num', how='left')
    all_df['engraftment_frequency'] = all_df['engraftment_frequency'].fillna(0)
    all_df.drop(columns=['sgb_num'], inplace=True)

    # Save full merged output
    full_output_csv = output_csv.replace(".csv", "_with_nonengrafted.csv")
    all_df.to_csv(full_output_csv, index=False)
    print(f"Full data (including non-engrafted SGBs) saved to {full_output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Associate engraftment counts with MetaPhlAn relative abundances")
    parser.add_argument('--engraftment_csv', required=True, help="CSV file with engraftment frequencies")
    parser.add_argument('--metaphlan_tsv', required=True, help="MetaPhlAn TSV file")
    parser.add_argument('--output_csv', required=True, help="Output CSV filename for merged data (engrafted only)")
    parser.add_argument('--sample_list', required=True, help="Text file with sample names to include (one per line)")

    args = parser.parse_args()
    main(args.engraftment_csv, args.metaphlan_tsv, args.sample_list, args.output_csv)
