import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import argparse
import numpy as np


def plot_relationship(data_csv, output_prefix):
    # Load merged data
    df = pd.read_csv(data_csv)

    # Basic checks
    if df.empty:
        print("Input data is empty!")
        return

    # Extract variables
    x = df['relative_abundance']
    y = df['engraftment_fraction_global']

    # Calculate Spearman correlation
    r, p = spearmanr(x, y)
    print(f"Spearman correlation: r = {r:.3f}, p = {p:.3g}")

    # Scatter plot with regression line (using LOWESS smoothing)
    plt.figure(figsize=(8,6))
    sns.regplot(x=x, y=y, scatter_kws={'s':30}, lowess=True, line_kws={'color':'red'})
    plt.xlabel('Relative Abundance')
    plt.ylabel('Engraftment Fraction (Global)')
    plt.title(f'Spearman r={r:.3f}, p={p:.3g}')
    plt.savefig(f'{output_prefix}_scatter.png')
    plt.close()

    # Histograms of each variable
    plt.figure(figsize=(12,5))

    plt.subplot(1,2,1)
    sns.histplot(x, kde=True, bins=30, color='blue')
    plt.xlabel('Relative Abundance')
    plt.title('Distribution of Relative Abundance')

    plt.subplot(1,2,2)
    sns.histplot(y, kde=True, bins=30, color='green')
    plt.xlabel('Engraftment Fraction (Global)')
    plt.title('Distribution of Engraftment Fraction')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_histograms.png')
    plt.close()

    # Optional: log-log scatter plot (if data has many zeros, handle carefully)
    # Add a small pseudocount to avoid log(0)
    x_log = x + 1e-6
    y_log = y + 1e-6

    plt.figure(figsize=(8,6))
    sns.scatterplot(x=np.log10(x_log), y=np.log10(y_log), s=30)
    plt.xlabel('Log10(Relative Abundance + 1e-6)')
    plt.ylabel('Log10(Engraftment Fraction + 1e-6)')
    plt.title('Log-Log Scatter Plot')
    plt.savefig(f'{output_prefix}_loglog_scatter.png')
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Spearman correlation and distributions between relative abundance and engraftment.")
    parser.add_argument('--merged_csv', required=True, help="CSV with merged data of relative abundance and engraftment_fraction_global")
    parser.add_argument('--output_prefix', required=True, help="Prefix for output plot files")
    args = parser.parse_args()

    plot_relationship(args.merged_csv, args.output_prefix)
