import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
import argparse
import numpy as np
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

# Set style for better-looking plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def plot_relationship(data_csv, output_prefix):
    # Load merged data
    df = pd.read_csv(data_csv)
    
    # Basic checks
    if df.empty:
        print("❌ Input data is empty!")
        return
    
    print(f"Loaded data: {len(df)} SGBs")
    
    # Extract variables
    x = df['relative_abundance']
    y = df['engraftment_fraction_global']
    
    # Remove any rows with NaN values
    valid_mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[valid_mask]
    y_clean = y[valid_mask]
    
    if len(x_clean) == 0:
        print("❌ No valid data points after removing NaN values!")
        return
    
    print(f"Valid data points: {len(x_clean)}")
    
    # Calculate correlations
    spearman_r, spearman_p = spearmanr(x_clean, y_clean)
    pearson_r, pearson_p = pearsonr(x_clean, y_clean)
    
    print(f"Spearman correlation: r = {spearman_r:.3f}, p = {spearman_p:.3e}")
    print(f"Pearson correlation: r = {pearson_r:.3f}, p = {pearson_p:.3e}")
    
    # Summary statistics
    print(f"Data summary:")
    print(f"   Relative abundance: {x_clean.min():.6f} - {x_clean.max():.6f} (median: {x_clean.median():.6f})")
    print(f"   Engraftment fraction: {y_clean.min():.3f} - {y_clean.max():.3f} (median: {y_clean.median():.3f})")
    
    # Create a comprehensive figure with multiple subplots
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Main scatter plot with regression line
    ax1 = plt.subplot(2, 3, (1, 2))
    
    # Color points by density for better visualization
    from scipy.stats import gaussian_kde
    if len(x_clean) > 10:
        xy = np.vstack([x_clean, y_clean])
        density = gaussian_kde(xy)(xy)
        scatter = ax1.scatter(x_clean, y_clean, c=density, s=40, alpha=0.7, cmap='viridis')
        plt.colorbar(scatter, ax=ax1, label='Point Density')
    else:
        ax1.scatter(x_clean, y_clean, s=40, alpha=0.7)
    
    # Add regression line
    if len(x_clean) > 3:
        sns.regplot(x=x_clean, y=y_clean, scatter=False, lowess=True, 
                   line_kws={'color':'red', 'linewidth': 2}, ax=ax1)
    
    ax1.set_xlabel('Relative Abundance', fontsize=12)
    ax1.set_ylabel('Engraftment Fraction', fontsize=12)
    ax1.set_title(f'SGB Abundance vs Engraftment\nSpearman r={spearman_r:.3f} (p={spearman_p:.2e})', 
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # 2. Log-log scatter plot (improved)
    ax2 = plt.subplot(2, 3, 3)
    
    # Handle zeros more intelligently
    x_nonzero = x_clean[x_clean > 0]
    y_nonzero = y_clean[x_clean > 0]
    y_nonzero = y_nonzero[y_nonzero > 0]
    x_nonzero = x_nonzero[y_clean[x_clean > 0] > 0]
    
    if len(x_nonzero) > 0:
        ax2.scatter(np.log10(x_nonzero), np.log10(y_nonzero), s=30, alpha=0.7)
        if len(x_nonzero) > 3:
            sns.regplot(x=np.log10(x_nonzero), y=np.log10(y_nonzero), scatter=False, 
                       line_kws={'color':'red'}, ax=ax2)
        ax2.set_xlabel('Log₁₀(Relative Abundance)', fontsize=10)
        ax2.set_ylabel('Log₁₀(Engraftment Fraction)', fontsize=10)
        ax2.set_title(f'Log-Log Plot\n({len(x_nonzero)} non-zero points)', fontsize=12)
    else:
        ax2.text(0.5, 0.5, 'No non-zero\ndata points', ha='center', va='center', 
                transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Log-Log Plot', fontsize=12)
    
    ax2.grid(True, alpha=0.3)
    
    # 3. Relative abundance distribution
    ax3 = plt.subplot(2, 3, 4)
    
    # Use log scale if data spans multiple orders of magnitude
    if x_clean.max() / x_clean[x_clean > 0].min() > 100:
        bins = np.logspace(np.log10(x_clean[x_clean > 0].min()), 
                          np.log10(x_clean.max()), 30)
        ax3.hist(x_clean[x_clean > 0], bins=bins, alpha=0.7, color='skyblue', edgecolor='black')
        ax3.set_xscale('log')
        ax3.set_xlabel('Relative Abundance (log scale)', fontsize=10)
    else:
        ax3.hist(x_clean, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
        ax3.set_xlabel('Relative Abundance', fontsize=10)
    
    ax3.set_ylabel('Frequency', fontsize=10)
    ax3.set_title(f'Abundance Distribution\n(n={len(x_clean)})', fontsize=12)
    ax3.grid(True, alpha=0.3)
    
    # 4. Engraftment distribution
    ax4 = plt.subplot(2, 3, 5)
    ax4.hist(y_clean, bins=30, alpha=0.7, color='lightcoral', edgecolor='black')
    ax4.set_xlabel('Engraftment Fraction', fontsize=10)
    ax4.set_ylabel('Frequency', fontsize=10)
    ax4.set_title(f'Engraftment Distribution\n(n={len(y_clean)})', fontsize=12)
    ax4.grid(True, alpha=0.3)
    
    # 5. Summary statistics table
    ax5 = plt.subplot(2, 3, 6)
    ax5.axis('off')
    
    # Create summary table
    stats_data = [
        ['Metric', 'Rel. Abundance', 'Engraftment Frac.'],
        ['Count', f'{len(x_clean)}', f'{len(y_clean)}'],
        ['Mean', f'{x_clean.mean():.4f}', f'{y_clean.mean():.3f}'],
        ['Median', f'{x_clean.median():.4f}', f'{y_clean.median():.3f}'],
        ['Std Dev', f'{x_clean.std():.4f}', f'{y_clean.std():.3f}'],
        ['Min', f'{x_clean.min():.4f}', f'{y_clean.min():.3f}'],
        ['Max', f'{x_clean.max():.4f}', f'{y_clean.max():.3f}'],
        ['', '', ''],
        ['Correlations', 'r-value', 'p-value'],
        ['Spearman', f'{spearman_r:.3f}', f'{spearman_p:.2e}'],
        ['Pearson', f'{pearson_r:.3f}', f'{pearson_p:.2e}']
    ]
    
    # Create table
    table = ax5.table(cellText=stats_data[1:], colLabels=stats_data[0],
                     cellLoc='center', loc='center',
                     colWidths=[0.3, 0.35, 0.35])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)
    
    # Style the table
    for i in range(len(stats_data)):
        for j in range(len(stats_data[0])):
            cell = table[(i, j)]
            if i == 0:  # Header
                cell.set_facecolor('#4CAF50')
                cell.set_text_props(weight='bold', color='white')
            elif i == 8:  # Correlation header
                cell.set_facecolor('#2196F3')
                cell.set_text_props(weight='bold', color='white')
            elif i in [9, 10]:  # Correlation values
                cell.set_facecolor('#E3F2FD')
    
    ax5.set_title('Summary Statistics', fontsize=12, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_comprehensive_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a focused scatter plot for publication
    plt.figure(figsize=(8, 6))
    
    # High-quality scatter plot
    if len(x_clean) > 10:
        xy = np.vstack([x_clean, y_clean])
        density = gaussian_kde(xy)(xy)
        scatter = plt.scatter(x_clean, y_clean, c=density, s=50, alpha=0.8, 
                            cmap='viridis', edgecolors='white', linewidth=0.5)
        plt.colorbar(scatter, label='Point Density')
    else:
        plt.scatter(x_clean, y_clean, s=50, alpha=0.8, edgecolors='white', linewidth=0.5)
    
    # Add regression line with confidence interval
    if len(x_clean) > 3:
        sns.regplot(x=x_clean, y=y_clean, scatter=False, lowess=True,
                   line_kws={'color':'red', 'linewidth': 2.5})
    
    plt.xlabel('Relative Abundance', fontsize=14, fontweight='bold')
    plt.ylabel('Engraftment Fraction', fontsize=14, fontweight='bold')
    plt.title(f'SGB Abundance vs Engraftment Efficiency\nSpearman ρ = {spearman_r:.3f}, p = {spearman_p:.2e}', 
              fontsize=16, fontweight='bold', pad=20)
    
    # Add text box with key statistics
    textstr = f'n = {len(x_clean)}\nSpearman ρ = {spearman_r:.3f}\np-value = {spearman_p:.2e}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
             verticalalignment='top', bbox=props)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_publication_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Generate summary report
    with open(f'{output_prefix}_analysis_report.txt', 'w') as f:
        f.write("SGB Abundance vs Engraftment Analysis Report\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Dataset: {data_csv}\n")
        f.write(f"Total SGBs analyzed: {len(x_clean)}\n\n")
        
        f.write("CORRELATION ANALYSIS:\n")
        f.write(f"Spearman correlation: ρ = {spearman_r:.4f}, p = {spearman_p:.4e}\n")
        f.write(f"Pearson correlation: r = {pearson_r:.4f}, p = {pearson_p:.4e}\n\n")
        
        f.write("RELATIVE ABUNDANCE STATISTICS:\n")
        f.write(f"Mean: {x_clean.mean():.6f}\n")
        f.write(f"Median: {x_clean.median():.6f}\n")
        f.write(f"Standard deviation: {x_clean.std():.6f}\n")
        f.write(f"Range: {x_clean.min():.6f} - {x_clean.max():.6f}\n\n")
        
        f.write("ENGRAFTMENT FRACTION STATISTICS:\n")
        f.write(f"Mean: {y_clean.mean():.4f}\n")
        f.write(f"Median: {y_clean.median():.4f}\n")
        f.write(f"Standard deviation: {y_clean.std():.4f}\n")
        f.write(f"Range: {y_clean.min():.4f} - {y_clean.max():.4f}\n\n")
        
        # Interpretation
        f.write("INTERPRETATION:\n")
        if abs(spearman_r) > 0.7:
            strength = "strong"
        elif abs(spearman_r) > 0.5:
            strength = "moderate"
        elif abs(spearman_r) > 0.3:
            strength = "weak"
        else:
            strength = "very weak"
        
        direction = "positive" if spearman_r > 0 else "negative"
        significance = "significant" if spearman_p < 0.05 else "not significant"
        
        f.write(f"There is a {strength} {direction} correlation between SGB relative abundance\n")
        f.write(f"and engraftment fraction that is {significance} (p = {spearman_p:.4e}).\n")
    
    print(f"\n📁 Generated files:")
    print(f"   • {output_prefix}_comprehensive_analysis.png - Detailed multi-panel analysis")
    print(f"   • {output_prefix}_publication_scatter.png - High-quality scatter plot")
    print(f"   • {output_prefix}_analysis_report.txt - Summary report")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Enhanced visualization of SGB abundance vs engraftment relationship")
    parser.add_argument('--merged_csv', required=True, help="CSV with merged data of relative abundance and engraftment_fraction_global")
    parser.add_argument('--output_prefix', required=True, help="Prefix for output plot files")
    args = parser.parse_args()

    plot_relationship(args.merged_csv, args.output_prefix)
