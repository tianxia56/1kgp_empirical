import pandas as pd
import os
from functools import reduce
import resource
import gc

# Define paths and populations
base_path = "/home/tx56/palmer_scratch/deepsweep_empirical"
daf_path = os.path.join(base_path, "1kgp_daf")
fst_path = os.path.join(base_path, "1kgp_fst")
normed_path = os.path.join(base_path, "normed")
output_path = os.path.join(base_path, "outputs")
target_pops = ["YRI", "CEU", "CHB", "BEB"]
pop_pairs = [("YRI", "CEU"), ("YRI", "CHB"), ("YRI", "BEB"), ("CEU", "CHB"), ("CEU", "BEB"), ("CHB", "BEB")]

# Create output directory if not exists
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Function to print memory usage
def print_memory_usage(stage):
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"Memory usage at {stage}: {mem} MB")

# Clear DataFrames to release memory
def clear_dfs(*dfs):
    for df in dfs:
        del df
    gc.collect()

# Process a single chromosome for each population
def process_chromosome(chr_num):
    for pop in target_pops:
        daf_file = os.path.join(daf_path, f"{pop}.{chr_num}.AA.tsv")
        empirical_file = os.path.join(output_path, f"{pop}.{chr_num}.empirical.tsv")
        if not os.path.exists(daf_file):
            continue

        df_empirical = pd.read_csv(daf_file, sep='\t').astype({'pos': 'Int64'})

        # Process FST files
        df_fst_combined = []
        for pair in pop_pairs:
            if pop in pair:
                fst_file = os.path.join(fst_path, f"{pair[0]}.{pair[1]}.ALL.chr{chr_num}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.weir.fst")
                if os.path.exists(fst_file):
                    df_fst = pd.read_csv(fst_file, sep='\t')
                    df_fst = df_fst[df_fst['WEIR_AND_COCKERHAM_FST'].notna()]
                    df_fst_combined.append(df_fst)

        if len(df_fst_combined) >= 3:
            df_fst_merged = reduce(lambda left, right: pd.merge(left, right, on='POS', how='outer'), df_fst_combined)
            fst_columns = [col for col in df_fst_merged.columns if col.startswith('WEIR_AND_COCKERHAM_FST')]
            df_fst_merged['mean_fst'] = df_fst_merged[fst_columns].mean(axis=1)
            df_fst_merged['mean_fst'] = df_fst_merged['mean_fst'].apply(lambda x: round(x, 5))
            df_fst_final = df_fst_merged[['POS', 'mean_fst']].drop_duplicates()
            df_empirical = pd.merge(df_empirical, df_fst_final.rename(columns={"POS": "pos"}), on='pos', how='outer')
            clear_dfs(df_fst_combined, df_fst_merged, df_fst_final)

        # Process NSL files
        nsl_file = os.path.join(normed_path, f"norm_nsl_{pop}_{chr_num}.tsv")
        if os.path.exists(nsl_file):
            df_nsl = pd.read_csv(nsl_file, sep='\t').drop_duplicates(subset='pos').rename(columns={'rs': 'rsid_nsl'}).astype({'pos': 'Int64'})
            df_empirical = pd.merge(df_empirical, df_nsl[['pos', 'rsid_nsl', 'norm_nsl']], on='pos', how='outer')
            clear_dfs(df_nsl)

        # Process IHH12 files
        ihh12_file = os.path.join(normed_path, f"norm_ihh12_{pop}_{chr_num}.tsv")
        if os.path.exists(ihh12_file):
            df_ihh12 = pd.read_csv(ihh12_file, sep='\t').drop_duplicates(subset='pos').rename(columns={'id': 'rsid_ihh12'}).astype({'pos': 'Int64'})
            df_empirical = pd.merge(df_empirical, df_ihh12[['pos', 'rsid_ihh12', 'norm_ihh12']], on='pos', how='outer')
            clear_dfs(df_ihh12)

        # Create inclusive rsid column
        df_empirical['rsid_inclusive'] = df_empirical[['rsid_nsl', 'rsid_ihh12']].bfill(axis=1).iloc[:, 0]

        # Process DELIHH files
        delihh_file = os.path.join(normed_path, f"norm_delihh_{pop}_{chr_num}.tsv")
        if os.path.exists(delihh_file):
            df_delihh = pd.read_csv(delihh_file, sep='\t').drop_duplicates(subset='Index').rename(columns={'ID': 'rsid_delihh'}).astype({'Index': 'Int64'})
            df_empirical = pd.merge(df_empirical, df_delihh[['rsid_delihh', 'norm_iHS', 'norm_delihh']], left_on='rsid_inclusive', right_on='rsid_delihh', how='outer')
            df_empirical.drop(columns=['rsid_delihh'], inplace=True, errors='ignore')
            clear_dfs(df_delihh)
        
        df_empirical.to_csv(empirical_file, sep='\t', index=False)

        # Process XPEHH files
        norm_xpehh_files = [
            os.path.join(normed_path, f)
            for f in os.listdir(normed_path)
            if f.startswith(f"norm_XPEHH_{pop}_") and f"_{chr_num}.tsv" in f
        ]
        
        if norm_xpehh_files:
            df_xpehh_combined = []
            for i, xpehh_file in enumerate(norm_xpehh_files):
                if os.path.exists(xpehh_file):
                    df_xpehh = pd.read_csv(xpehh_file, sep='\t')
                    df_xpehh = df_xpehh.drop_duplicates(subset='Index').rename(columns={'ID': 'rsid', 'norm_XPEHH': f'norm_XPEHH_{i}'}).astype({'Index': 'Int64'})
                    df_xpehh_combined.append(df_xpehh)

            # Use inner join for calculating max_xpehh
            df_xpehh_merged = reduce(lambda left, right: pd.merge(left, right, on='rsid', how='inner'), df_xpehh_combined)
            xpehh_columns = [col for col in df_xpehh_merged.columns if col.startswith('norm_XPEHH_')]
            df_xpehh_merged['max_xpehh'] = df_xpehh_merged[xpehh_columns].apply(lambda row: row.max(), axis=1)

            df_empirical = pd.merge(df_empirical, df_xpehh_merged[['Index', 'rsid', 'max_xpehh']], left_on='rsid_inclusive', right_on='rsid', how='outer')
            df_empirical.drop(columns=['Index', 'rsid', 'rsid_nsl', 'rsid_ihh12'], inplace=True, errors='ignore')
            df_empirical.rename(columns={'rsid_inclusive': 'rsid'}, inplace=True)
            df_empirical.to_csv(empirical_file, sep='\t', index=False)
            clear_dfs(df_xpehh_combined, df_xpehh_merged)

        # Convert mean_fst to 5 significant figures, replace NaNs with pd.NA, reorder rows and columns
        df_empirical['mean_fst'] = df_empirical['mean_fst'].apply(lambda x: round(x, 5))
        df_empirical = df_empirical.sort_values(by='pos')
        df_empirical = df_empirical.fillna(pd.NA)
        df_empirical = df_empirical[['rsid'] + [col for col in df_empirical.columns if col != 'rsid']]
        
        # Save final DataFrame
        df_empirical.to_csv(empirical_file, sep='\t', index=False)
        clear_dfs(df_empirical)
        gc.collect()

# Process each chromosome sequentially
for chr_num in range(1, 23):
    process_chromosome(chr_num)

print("Processing completed.")
