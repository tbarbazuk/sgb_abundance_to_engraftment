import pandas as pd
from collections import Counter
import argparse
import sys

def filter_by_timepoint(metadata_csv, timepoint='t7'):
    df = pd.read_csv(metadata_csv, sep='\t')
    filtered_df = df[df['timepoint'] == timepoint].copy()
    print(f"Found {len(filtered_df)} samples at timepoint '{timepoint}'")
    return filtered_df

def get_sample_to_donor_map(filtered_df):
    sorted_df = filtered_df.sort_values(by='donor_sample_name')
    total = len(sorted_df)
    nan_count = sorted_df['donor_sample_name'].isna().sum()
    fraction_missing = nan_count / total if total > 0 else 0
    sample_to_donor = dict(zip(sorted_df['sample_id'], sorted_df['donor_sample_name']))
    print(f"Number of samples missing donor_sample_name: {nan_count}")
    print(f"Fraction missing: {fraction_missing:.2%}")
    return sample_to_donor, nan_count, fraction_missing

def count_samples_per_donor(sample_to_donor):
    """
    Counts the number of samples that received an FMT from each donor.

    Parameters:
    - sample_to_donor (dict): A dictionary mapping sample_id -> donor_sample_name.

    Returns:
    - donor_counts (dict): A dictionary mapping donor_sample_name -> count of samples.
    """
    donor_counts = Counter(sample_to_donor.values())
    
    print("Number of samples per donor:")
    for donor, count in donor_counts.items():
        print(f"  {donor}: {count}")
    
    return donor_counts

def safe_strain_list(val):
    if pd.isna(val) or val.strip() == "[]":
        return set()
    val = val.strip("[]").replace("'", "").replace('"', "")
    return set([s.strip() for s in val.split() if s.strip()])

def get_matched_samples(filtered_df):
    matched_df = filtered_df[filtered_df['donor_sample_name'].notna()].copy()
    matched_df = matched_df.sort_values(by='donor_sample_name')
    matched_map = dict(zip(matched_df['sample_id'], matched_df['donor_sample_name']))
    print(f"\nMatched samples (excluding placebo): {len(matched_df)} found.")
    return matched_df, matched_map

def load_engraftment_data(filepath, timepoint_filter='t7.0'):
    df = pd.read_csv(filepath, index_col=0)
    print(f"Engraftment data shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    filtered_df = df[df['timepoint'] == timepoint_filter].copy()
    result = filtered_df[['subject_id', 'eng_strains', 'total_donor_strains']].copy()
    result.index.name = 'sample_id'
    print(f"Filtered engraftment rows for {timepoint_filter}: {result.shape[0]}")
    return result

def find_strain_matches(df, eng_col='eng_strains', donor_col='total_donor_strains'):
    matches = []
    for _, row in df.iterrows():
        eng_strains = safe_strain_list(row[eng_col])
        donor_strains = safe_strain_list(row[donor_col])
        matched = eng_strains.intersection(donor_strains)
        matches.append(list(sorted(matched)) if matched else [])
    df['matched_strains'] = matches
    return df

def count_engrafted_sgbs(df, match_col='matched_strains'):
    df['matched_strains_list'] = df[match_col].apply(lambda x: x if isinstance(x, list) else [])
    exploded_df = df.explode('matched_strains_list')
    exploded_df = exploded_df[exploded_df['matched_strains_list'].notna() & (exploded_df['matched_strains_list'] != '')]
    sgb_counts = (
        exploded_df
        .groupby(['donor_sample_name', 'matched_strains_list'])
        .size()
        .reset_index(name='count')
        .rename(columns={'matched_strains_list': 'sgb_id'})
    )
    return sgb_counts

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Analyze engraftment data for FMT samples at specified timepoints',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python script.py t7 metadata.tsv engraftment.csv
  python script.py t30 metadata_ss_pipeline_TACITO_fixed.tsv eng.csv
        """
    )
    
    parser.add_argument(
        'timepoint',
        help='Timepoint to analyze (e.g., t7, t30)'
    )
    
    parser.add_argument(
        'metadata_file',
        help='Path to the metadata TSV file'
    )
    
    parser.add_argument(
        'engraftment_file',
        help='Path to the engraftment CSV file'
    )
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    timepoint_of_interest = args.timepoint
    metadata_path = args.metadata_file
    engraftment_path = args.engraftment_file
    
    # Convert timepoint for engraftment data (add .0 if not present)
    engraftment_timepoint = f"{timepoint_of_interest}.0" if not timepoint_of_interest.endswith('.0') else timepoint_of_interest
    
    print(f"Analyzing timepoint: {timepoint_of_interest}")
    print(f"Metadata file: {metadata_path}")
    print(f"Engraftment file: {engraftment_path}")
    print(f"Engraftment timepoint filter: {engraftment_timepoint}")
    print("-" * 50)

    try:
        # Step 1: Filter metadata
        filtered_metadata = filter_by_timepoint(metadata_path, timepoint=timepoint_of_interest)

        # Step 2: Create mapping and get stats
        sample_to_donor_map, missing_count, missing_fraction = get_sample_to_donor_map(filtered_metadata)
        
        # Step 3: Get matched samples
        matched_df, matched_map = get_matched_samples(filtered_metadata)

        # Step 4: Load engraftment data
        engraftment_df = load_engraftment_data(engraftment_path, timepoint_filter=engraftment_timepoint)

        # Step 5: Add donor_sample_name to engraftment data
        matched_df['subject_id'] = matched_df['sample_id'].str.extract(r'FMT_(T\d{3})')
        subject_to_donor_df = matched_df[['subject_id', 'donor_sample_name']].drop_duplicates()
        engraftment_df = engraftment_df.reset_index()
        engraftment_df = pd.merge(engraftment_df, subject_to_donor_df, on='subject_id', how='left')
        # Step 5 continued: Remove donors starting with 'DT005'
        engraftment_df = engraftment_df[~engraftment_df['donor_sample_name'].str.startswith("DTT005", na=False)]

        unmatched = engraftment_df['donor_sample_name'].isna().sum()
        print(f"\nDonor assignment complete. Unmatched entries: {unmatched}")

        # Step 6: Match strains from donor and recipient
        engraftment_df = find_strain_matches(engraftment_df)

        # Step 7: Count engrafted SGBs per donor
        sgb_count_df = count_engrafted_sgbs(engraftment_df)
        
        # Count total recipient samples per donor at this timepoint
        samples_per_donor = (
            engraftment_df
            .groupby('donor_sample_name')['sample_id']
            .nunique()
            .reset_index(name='total_samples')
        )

        print("\nSamples per donor (used as denominator for per-donor engraftment fraction):")
        print(samples_per_donor)

        # Merge with the SGB counts
        sgb_fraction_df = pd.merge(sgb_count_df, samples_per_donor, on='donor_sample_name', how='left')

        # Compute the per-donor engraftment fraction
        sgb_fraction_df['engraftment_fraction'] = sgb_fraction_df['count'] / sgb_fraction_df['total_samples']
        
        # --- Compute GLOBAL denominator engraftment fraction ---
        total_fmt_samples = engraftment_df['donor_sample_name'].notna().sum()
        print(f"\nTotal FMT samples at {timepoint_of_interest} (excluding placebo): {total_fmt_samples}")

        sgb_fraction_global_df = sgb_count_df.copy()
        sgb_fraction_global_df['engraftment_fraction_global'] = (
            sgb_fraction_global_df['count'] / total_fmt_samples
        )
        
        # Export global engraftment fractions CSV:
        output_filename = f"engraftment_frequencies_{timepoint_of_interest}.csv"
        sgb_fraction_global_df.to_csv(output_filename, index=False)
        print(f"\nEngraftment frequency table with global denominator exported to {output_filename}")

        # Get top 10 SGBs per donor using global denominator
        top_global_fraction_sgbs = (
            sgb_fraction_global_df
            .sort_values(['donor_sample_name', 'engraftment_fraction_global'], ascending=[True, False])
            .groupby('donor_sample_name')
            .head(10)
        )

        # Sort for readability
        sgb_fraction_df = sgb_fraction_df.sort_values(['donor_sample_name', 'engraftment_fraction'], ascending=[True, False])

        # Step 8: Display top 10 SGBs per donor
        top_sgbs = (
            sgb_count_df
            .sort_values(['donor_sample_name', 'count'], ascending=[True, False])
            .groupby('donor_sample_name')
            .head(10)
        )
        
        # Step 9: Display top 10 SGBs per donor by engraftment fraction
        top_fraction_sgbs = (
            sgb_fraction_df
            .sort_values(['donor_sample_name', 'engraftment_fraction'], ascending=[True, False])  # <-- required
            .groupby('donor_sample_name')
            .head(10)
        )
        
        print("Sample to donor map")
        print(sample_to_donor_map)
        
        print("\n--- Total recipient samples per donor (from sample_to_donor_map) ---")
        print(samples_per_donor.sort_values('total_samples', ascending=False))
        
        print("\n--- Engraftment DataFrame (after donor matching) ---")
        print(engraftment_df[['sample_id', 'subject_id', 'donor_sample_name', 'eng_strains', 'total_donor_strains']].head(5))

        print("\n--- Full SGB count DataFrame ---")
        print(sgb_count_df.head(10)) 

        print("\n--- Raw types and contents of eng_strains and total_donor_strains ---")
        print(engraftment_df[['eng_strains', 'total_donor_strains']].head(3))
        print("\nTypes:")
        print(type(engraftment_df['eng_strains'].iloc[0]), type(engraftment_df['total_donor_strains'].iloc[0]))

        print("\nTop SGBs by engraftment fraction per donor:")
        print(top_fraction_sgbs[['donor_sample_name', 'sgb_id', 'count', 'total_samples', 'engraftment_fraction']])

        print("\nSGBs by engraftment fraction across all donor samples (GLOBAL denominator):")
        print(top_global_fraction_sgbs[['donor_sample_name', 'sgb_id', 'count', 'engraftment_fraction_global']])

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

# === MAIN EXECUTION ===
if __name__ == "__main__":
    main()
