import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Define the positions and rsids of the selected variants
selected_variants = {
    "LCT": {"pos": 136608646, "rsid": "rs4988235"},
    "SLC24A5": {"pos": 48533135, "rsid": "rs1426654"},
    "EDAR": [
        {"pos": 108348155, "rsid": "rs3827760"},
        {"pos": 108418833, "rsid": "rs146567337"}
    ],
    "DARC": {"pos": 159174698, "rsid": "rs2814778"},
    "SLC45A2": {"pos": 33951691, "rsid": "rs16891982"}  # Added SLC45A2
}

# Define the coding regions of the genes
coding_regions = {
    "LCT": (136545000, 136610000),
    "SLC24A5": (48520000, 48540000),
    "EDAR": (108340000, 108420000),
    "DARC": (159160000, 159180000),
    "SLC45A2": (33950000, 33960000)  # Added coding region for SLC45A2
}

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
    columns_to_plot = ['Fst', 'ΔDAF', 'iSAFE', 'iHS', 'nSL', 'iHH12', 'ΔiHH', 'XPEHH']

    # Check if all columns to plot exist in the data
    missing_columns = [col for col in columns_to_plot if col not in data.columns]
    if missing_columns:
        print(f"Missing columns in data: {missing_columns}")
        return

    # Extract gene and population from the file name
    file_name = os.path.basename(file_path)
    gene, pop = file_name.split('.')[0], file_name.split('.')[1]

    # Create a combined plot with shared x-axis and 16:9 aspect ratio
    fig, axes = plt.subplots(len(columns_to_plot) + 1, 1, figsize=(16, 12), sharex=True)

    for i, column in enumerate(columns_to_plot):
        sns.scatterplot(x='pos', y=column, data=data, ax=axes[i], color='black', s=10, edgecolor=None)
        axes[i].set_ylabel(column, rotation=0, labelpad=40, fontsize=16)
        axes[i].tick_params(axis='y', which='both', length=0)  # Remove y-axis ticks
        axes[i].spines['top'].set_visible(False)
        axes[i].spines['right'].set_visible(False)
        axes[i].spines['left'].set_visible(False)

        # Add vertical red dashed lines at the selected variant positions
        if gene in selected_variants:
            if isinstance(selected_variants[gene], list):
                for variant in selected_variants[gene]:
                    axes[i].axvline(x=variant["pos"], color='red', linestyle='--')
            else:
                selected_pos = selected_variants[gene]["pos"]
                axes[i].axvline(x=selected_pos, color='red', linestyle='--')

    # Add the coding region bar plot
    if gene in coding_regions:
        coding_start, coding_end = coding_regions[gene]
        axes[-1].barh(y=0, width=coding_end - coding_start, left=coding_start, height=0.5, color='blue')
        axes[-1].set_yticks([])
        axes[-1].set_ylabel('Coding Region', rotation=0, labelpad=40, fontsize=16)
        axes[-1].spines['top'].set_visible(False)
        axes[-1].spines['right'].set_visible(False)
        axes[-1].spines['left'].set_visible(False)
        axes[-1].spines['bottom'].set_visible(False)

    axes[-1].set_xlabel('Position', fontsize=18)
    axes[-1].ticklabel_format(style='plain', axis='x')  # Avoid scientific notation

    # Add the selected variant rsid and position to the main title
    if gene in selected_variants:
        if isinstance(selected_variants[gene], list):
            variant_info = ', '.join([f'{v["rsid"]} at {v["pos"]}' for v in selected_variants[gene]])
        else:
            variant_info = f'{selected_variants[gene]["rsid"]} at {selected_variants[gene]["pos"]}'
        plt.suptitle(f'{gene} - {pop} - {variant_info}', fontsize=20)
    else:
        plt.suptitle(f'{gene} - {pop}', fontsize=20)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=1)  # Adjust gaps between plots

    # Save the combined plot to a file
    output_plot_path = f"../sel_loci/{gene}.{pop}.png"
    plt.savefig(output_plot_path)
    plt.close()

    print(f"Combined plot saved to {output_plot_path}")

if __name__ == "__main__":
    stats_files = [
        "../sel_loci/LCT.CEU.stats",
        "../sel_loci/SLC24A5.BEB.stats",  # Changed to BEB
        "../sel_loci/EDAR.CHB.stats",
        "../sel_loci/DARC.YRI.stats",
        "../sel_loci/SLC45A2.CEU.stats"  # Added SLC45A2 for CEU
    ]
    
    for file_path in stats_files:
        plot_stats(file_path)
