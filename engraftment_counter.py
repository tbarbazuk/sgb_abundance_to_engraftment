import pandas as pd
from collections import Counter
import argparse
import sys
import os

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

def create_comprehensive_donor_csvs(engraftment_df, timepoint, prefix, output_dir):
    """
    Create comprehensive donor-specific CSV files with all engraftment information
    and frequency calculations.
    """
    print(f"\nCreating comprehensive donor-specific CSV files...")
    
    # Get unique donors
    donors = engraftment_df['donor_sample_name'].dropna().unique()
    
    for donor in sorted(donors):
        # Filter data for this donor
        donor_df = engraftment_df[engraftment_df['donor_sample_name'] == donor].copy()
        
        # Calculate basic stats for this donor
        total_samples_this_donor = len(donor_df)
        
        # Create detailed sample-level information
        sample_details = []
        
        for _, row in donor_df.iterrows():
            sample_info = {
                'sample_id': row['sample_id'],
                'subject_id': row['subject_id'],
                'donor_sample_name': row['donor_sample_name'],
                'timepoint': timepoint,
                'eng_strains': row['eng_strains'],
                'total_donor_strains': row['total_donor_strains'],
                'matched_strains': row['matched_strains'],
                'num_engrafted_strains': len(row['matched_strains']) if row['matched_strains'] else 0,
                'total_donor_strains_count': len(safe_strain_list(row['total_donor_strains'])),
                'engraftment_success': len(row['matched_strains']) > 0 if row['matched_strains'] else False
            }
            sample_details.append(sample_info)
        
        # Create sample-level DataFrame
        sample_df = pd.DataFrame(sample_details)
        
        # Calculate strain-level engraftment frequencies for this donor
        strain_stats = []
        
        # Get all strains that were engrafted for this donor
        all_engrafted_strains = []
        for matched_list in donor_df['matched_strains']:
            if matched_list:
                all_engrafted_strains.extend(matched_list)
        
        # Count occurrences of each strain
        strain_counts = Counter(all_engrafted_strains)
        
        for strain, count in strain_counts.items():
            strain_info = {
                'sgb_id': strain,
                'engraftment_count': count,
                'total_samples': total_samples_this_donor,
                'engraftment_frequency': count / total_samples_this_donor,
                'engraftment_percentage': (count / total_samples_this_donor) * 100
            }
            strain_stats.append(strain_info)
        
        # Create strain-level DataFrame
        strain_df = pd.DataFrame(strain_stats)
        strain_df = strain_df.sort_values('engraftment_frequency', ascending=False)
        
        # Calculate donor summary statistics
        summary_stats = {
            'donor_sample_name': donor,
            'timepoint': timepoint,
            'total_samples': total_samples_this_donor,
            'samples_with_engraftment': sample_df['engraftment_success'].sum(),
            'samples_without_engraftment': total_samples_this_donor - sample_df['engraftment_success'].sum(),
            'overall_engraftment_rate': sample_df['engraftment_success'].mean(),
            'mean_strains_per_sample': sample_df['num_engrafted_strains'].mean(),
            'median_strains_per_sample': sample_df['num_engrafted_strains'].median(),
            'max_strains_per_sample': sample_df['num_engrafted_strains'].max(),
            'min_strains_per_sample': sample_df['num_engrafted_strains'].min(),
            'total_unique_engrafted_strains': len(strain_counts),
            'most_frequent_strain': strain_df.iloc[0]['sgb_id'] if len(strain_df) > 0 else None,
            'highest_engraftment_frequency': strain_df.iloc[0]['engraftment_frequency'] if len(strain_df) > 0 else 0
        }
        
        # Create summary DataFrame
        summary_df = pd.DataFrame([summary_stats])
        
        # Save comprehensive CSV file
        donor_filename = f"{output_dir}/{prefix}{donor}_comprehensive_{timepoint}.csv"
        
        # Create a comprehensive output by combining all information
        with open(donor_filename, 'w') as f:
            # Write header
            f.write(f"# Comprehensive Engraftment Analysis for Donor: {donor}\n")
            f.write(f"# Timepoint: {timepoint}\n")
            f.write(f"# Generated by engraftment analysis script\n")
            f.write(f"# Total samples: {total_samples_this_donor}\n")
            f.write(f"# Samples with engraftment: {sample_df['engraftment_success'].sum()}\n")
            f.write(f"# Overall engraftment rate: {sample_df['engraftment_success'].mean():.2%}\n")
            f.write("\n")
            
            # Write summary statistics
            f.write("## DONOR SUMMARY STATISTICS\n")
            summary_df.to_csv(f, index=False)
            f.write("\n")
            
            # Write strain-level statistics
            f.write("## STRAIN-LEVEL ENGRAFTMENT FREQUENCIES\n")
            strain_df.to_csv(f, index=False)
            f.write("\n")
            
            # Write sample-level details
            f.write("## SAMPLE-LEVEL DETAILS\n")
            sample_df.to_csv(f, index=False)
        
        print(f"Saved comprehensive analysis for {donor} to {donor_filename}")
        
        # Also save separate files for easier analysis
        summary_file = f"{output_dir}/{prefix}{donor}_summary_{timepoint}.csv"
        strain_file = f"{output_dir}/{prefix}{donor}_strain_frequencies_{timepoint}.csv"
        sample_file = f"{output_dir}/{prefix}{donor}_sample_details_{timepoint}.csv"
        
        summary_df.to_csv(summary_file, index=False)
        strain_df.to_csv(strain_file, index=False)
        sample_df.to_csv(sample_file, index=False)
        
        print(f"  - Summary: {summary_file}")
        print(f"  - Strain frequencies: {strain_file}")
        print(f"  - Sample details: {sample_file}")

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

        # Step 8.1: Create comprehensive donor-specific CSV files (NEW IMPROVED VERSION)
        # Output directly to current directory instead of subdirectory
        create_comprehensive_donor_csvs(engraftment_df, timepoint_of_interest, prefix, ".")

        # Step 8.2: Write basic individual donor-specific CSVs (keeping original functionality)
        donor_dir = f"donor_engraftment_csvs_{timepoint_of_interest}"
        os.makedirs(donor_dir, exist_ok=True)

        for donor in sorted(remaining_donors):
            donor_df = engraftment_df[engraftment_df['donor_sample_name'] == donor].copy()
            donor_csv = f"{donor_dir}/{prefix}{donor}_engraftment.csv"
            donor_df.to_csv(donor_csv, index=False)
        
        print(f"\nBasic per-donor engraftment CSVs saved in folder: {donor_dir}")

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
        
        # Step 12.5: Generate per-donor engraftment CSVs (keeping original functionality)
        donor_sample_groups = engraftment_df.groupby('donor_sample_name')
        donor_csv_output_dir = f"donor_stats_{timepoint_of_interest}"
        os.makedirs(donor_csv_output_dir, exist_ok=True)

        for donor_name, group_df in donor_sample_groups:
            donor_sample_count = group_df.shape[0]
            donor_strains = (
                group_df.explode('matched_strains')
                .dropna(subset=['matched_strains'])
                .groupby('matched_strains')
                .size()
                .reset_index(name='count')
                .rename(columns={'matched_strains': 'sgb_id'})
            )
            donor_strains['engraftment_fraction'] = donor_strains['count'] / donor_sample_count

            donor_filename = f"{donor_csv_output_dir}/engraftment_{donor_name}_{timepoint_of_interest}.csv"
            donor_strains.to_csv(donor_filename, index=False)
            print(f"Exported per-donor engraftment stats to {donor_filename}")

        # Step 13: Print top 10 SGBs by global engraftment fraction
        top_global_fraction_sgbs = sgb_global_counts.sort_values(
            'engraftment_fraction_global', ascending=False
        ).head(10)

        print("\nTop 10 SGBs by global engraftment fraction:")
        print(top_global_fraction_sgbs)

        # Step 14: Print summary of all outputs
        print("\n" + "="*60)
        print("SUMMARY OF GENERATED FILES:")
        print("="*60)
        print(f"1. Global engraftment frequencies: {output_filename}")
        print(f"2. Donor list: {donor_output_file}")
        print(f"3. Basic donor CSVs: {donor_dir}/")
        print(f"4. Donor stats: {donor_csv_output_dir}/")
        print(f"5. COMPREHENSIVE donor analysis: {comprehensive_dir}/")
        print(f"   - Includes summary statistics, strain frequencies, and sample details")
        print(f"   - Each donor has 4 files: comprehensive, summary, strain frequencies, and sample details")

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

# === MAIN EXECUTION ===
if __name__ == "__main__":
    main()
