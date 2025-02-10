import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_stats(file_path):
    if not os.path.exists(file_path):
        print(f"Stats file {file_path} does not exist.")
        return

    # Load the data from the stats file
    data = pd.read_csv(file_path, sep='\t')

    # Rename the columns in the DataFrame
    data.rename(columns={
        'mean_fst': 'Fst',
        'deldaf': 'ΔDAF',
        'norm_iHS': 'iHS',  # Correct casing for norm_iHS
        'norm_nsl': 'nSL',
        'norm_ihh12': 'iHH12',
        'norm_delihh': 'ΔiHH',
        'max_xpehh': 'XPEHH'
    }, inplace=True)

    # Define the columns to plot (excluding 'maf' and 'daf')
    columns_to_plot = ['Fst', 'ΔDAF', 'iHS', 'nSL', 'iHH12', 'ΔiHH', 'XPEHH']

    # Check if all columns to plot exist in the data
    missing_columns = [col for col in columns_to_plot if col not in data.columns]
    if missing_columns:
        print(f"Missing columns in data: {missing_columns}")
        return

    # Extract gene and population from the file name
    file_name = os.path.basename(file_path)
    pop, chr = file_name.split('.')[0], file_name.split('.')[1]

    # Create a combined plot with shared x-axis and 16:9 aspect ratio
    fig, axes = plt.subplots(len(columns_to_plot), 1, figsize=(16, 12), sharex=True)

    for i, column in enumerate(columns_to_plot):
        sns.scatterplot(x='pos', y=column, data=data, ax=axes[i], color='black', s=10, edgecolor=None)
        axes[i].set_ylabel(column, rotation=0, labelpad=40, fontsize=16)
        axes[i].tick_params(axis='y', which='both', length=0)  # Remove y-axis ticks
        axes[i].spines['top'].set_visible(False)
        axes[i].spines['right'].set_visible(False)
        axes[i].spines['left'].set_visible(False)

    axes[-1].set_xlabel('Position', fontsize=18)
    axes[-1].ticklabel_format(style='plain', axis='x')  # Avoid scientific notation

    plt.suptitle(f'Chromosome {chr} - {pop}', fontsize=20)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=1)  # Adjust gaps between plots

    # Save the combined plot to a file
    output_plot_path = f"../sel_loci/{pop}.{chr}.png"
    plt.savefig(output_plot_path)
    plt.close()

    print(f"Combined plot saved to {output_plot_path}")

if __name__ == "__main__":
    stats_files = [
        "../outputs/BEB.22.empirical.tsv",
        "../outputs/CEU.22.empirical.tsv",
        "../outputs/CHB.22.empirical.tsv",
        "../outputs/YRI.22.empirical.tsv"
    ]
    
    for file_path in stats_files:
        plot_stats(file_path)
