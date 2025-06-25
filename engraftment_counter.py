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
    
    parser.add_argument(
        '--prefix',
        default='FMT_',
        help='Prefix to add to donor names in output file (default: FMT_)'
    )
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    timepoint_of_interest = args.timepoint
    metadata_path = args.metadata_file
    engraftment_path = args.engraftment_file
    prefix = args.prefix
    
    # Convert timepoint for engraftment data (add .0 if not present)
    engraftment_timepoint = f"{timepoint_of_interest}.0" if not timepoint_of_interest.endswith('.0') else timepoint_of_interest
    
    print(f"Analyzing timepoint: {timepoint_of_interest}")
    print(f"Metadata file: {metadata_path}")
    print(f"Engraftment file: {engraftment_path}")
    print(f"Engraftment timepoint filter: {engraftment_timepoint}")
    print(f"Prefix for donor names: '{prefix}'")
    print("-" * 50)

    try:
        # Step 1: Filter metadata
        filtered_metadata = filter_by_timepoint(metadata_path, timepoint=timepoint_of_interest)

        # Step 2: Create mapping and get stats
        sample_to_donor_map, missing_count, missing_fraction = get_sample_to_donor_map(filtered_metadata)
        print("sample_to_donor_map")
        print(sample_to_donor_map)        
        
        # Step 3: Get matched samples (samples with donor info)
        matched_df, matched_map = get_matched_samples(filtered_metadata)
        print(matched_map)
        
        # Step 4: Load engraftment data filtered by timepoint
        engraftment_df = load_engraftment_data(engraftment_path, timepoint_filter=engraftment_timepoint)

        # Step 5: Add donor_sample_name to engraftment data by merging on subject_id
        matched_df['subject_id'] = matched_df['sample_id'].str.extract(r'FMT_(T\d{3})')
        subject_to_donor_df = matched_df[['subject_id', 'donor_sample_name']].drop_duplicates()
        engraftment_df = engraftment_df.reset_index()
        engraftment_df = pd.merge(engraftment_df, subject_to_donor_df, on='subject_id', how='left')

        # Step 6: Remove donors starting with 'DTT005' (FILTERING STEP)
        print(f"\nBefore DTT005 filtering: {len(engraftment_df)} samples")
        engraftment_df = engraftment_df[~engraftment_df['donor_sample_name'].str.startswith("DTT005", na=False)]
        print(f"After DTT005 filtering: {len(engraftment_df)} samples")

        # Step 7: Generate donor text files AFTER filtering with prefix
        remaining_donors = engraftment_df['donor_sample_name'].dropna().unique()
        donor_output_file = f"donors_remaining_{timepoint_of_interest}.txt"
        
        print(f"\nFinal remaining donors after DTT005 filtering: {len(remaining_donors)}")
        print(f"Donor list (original): {sorted(remaining_donors)}")
        
        # Write donors with prefix to the text file
        with open(donor_output_file, 'w') as f:
            for donor in sorted(remaining_donors):
                prefixed_donor = f"{prefix}{donor}"
                f.write(f"{prefixed_donor}\n")
        
        # Show what was written to file
        prefixed_donors = [f"{prefix}{donor}" for donor in sorted(remaining_donors)]
        print(f"Donor list (with prefix '{prefix}'): {prefixed_donors}")
        print(f"Remaining donor sample names (with prefix) written to {donor_output_file}")

        unmatched = engraftment_df['donor_sample_name'].isna().sum()
        print(f"Donor assignment complete. Unmatched entries: {unmatched}")

        # Step 8: Find matched strains between donor and recipient
        engraftment_df = find_strain_matches(engraftment_df)

        # Step 9: Compute total FMT samples (exclude missing donors)
        total_fmt_samples = engraftment_df['donor_sample_name'].notna().sum()
        print(f"\nTotal FMT samples at {timepoint_of_interest} (excluding placebo and DTT005): {total_fmt_samples}")

        # Step 10: Count SGB occurrences globally (across all donors)
        sgb_global_counts = (
            engraftment_df.explode('matched_strains')
            .dropna(subset=['matched_strains'])
            .groupby('matched_strains')
            .size()
            .reset_index(name='count')
            .rename(columns={'matched_strains': 'sgb_id'})
        )

        # Step 11: Calculate global engraftment fraction
        sgb_global_counts['engraftment_fraction_global'] = sgb_global_counts['count'] / total_fmt_samples
        
        # Step 12: Export results
        output_filename = f"engraftment_frequencies_{timepoint_of_interest}.csv"
        sgb_global_counts.to_csv(output_filename, index=False)
        print(f"\nEngraftment frequency table exported to {output_filename}")

        # Step 13: Print top 10 SGBs by global engraftment fraction
        top_global_fraction_sgbs = sgb_global_counts.sort_values(
            'engraftment_fraction_global', ascending=False
        ).head(10)

        print("\nTop 10 SGBs by global engraftment fraction:")
        print(top_global_fraction_sgbs)

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

# === MAIN EXECUTION ===
if __name__ == "__main__":
    main()
