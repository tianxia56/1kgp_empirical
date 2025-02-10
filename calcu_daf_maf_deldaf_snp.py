import pandas as pd
import os
from functools import reduce
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

# Clear DataFrames to release memory
def clear_dfs(*dfs):
    for df in dfs:
        del df
    gc.collect()

# Function to preview the output TSV
def preview_output(df, step):
    print(f"Preview of output TSV after {step}:")
    print(df.head(2))

# Process a single chromosome for each population
def process_chromosome(chr_num):
    for pop in target_pops:
        daf_file = os.path.join(daf_path, f"{pop}.{chr_num}.AA.tsv")
        empirical_file = os.path.join(output_path, f"{pop}.{chr_num}.empirical.tsv")
        if not os.path.exists(daf_file):
            continue

        print(f"Reading DAF file for {pop} chromosome {chr_num}")
        df_empirical = pd.read_csv(daf_file, sep='\t').astype({'pos': 'Int64'})
        print(f"DAF file for {pop} chromosome {chr_num} read successfully")
        print(f"Word count for DAF file: {len(df_empirical)}")
        print(df_empirical.head(2))

        # Process FST files
        df_fst_combined = []
        for pair in pop_pairs:
            if pop in pair:
                fst_file = os.path.join(fst_path, f"{pair[0]}.{pair[1]}.ALL.chr{chr_num}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.weir.fst")
                if os.path.exists(fst_file):
                    print(f"Reading FST file for pair {pair} chromosome {chr_num}")
                    df_fst = pd.read_csv(fst_file, sep='\t')
                    df_fst = df_fst[df_fst['WEIR_AND_COCKERHAM_FST'].notna()]
                    df_fst_combined.append(df_fst)
                    print(f"FST file for pair {pair} chromosome {chr_num} read successfully")
                    print(f"Word count for FST file: {len(df_fst)}")
                    print(df_fst.head(2))

        if len(df_fst_combined) >= 3:
            print(f"Merging FST files for {pop} chromosome {chr_num}")
            df_fst_merged = reduce(lambda left, right: pd.merge(left, right, on='POS', how='outer'), df_fst_combined)
            print(f"Temporary merged FST file for {pop} chromosome {chr_num}")
            print(df_fst_merged.head(2))
            fst_columns = [col for col in df_fst_merged.columns if col.startswith('WEIR_AND_COCKERHAM_FST')]
            df_fst_merged['mean_fst'] = df_fst_merged[fst_columns].mean(axis=1)
            df_fst_merged['mean_fst'] = df_fst_merged['mean_fst'].apply(lambda x: round(x, 5))
            df_fst_final = df_fst_merged[['POS', 'mean_fst']].drop_duplicates()
            df_empirical = pd.merge(df_empirical, df_fst_final.rename(columns={"POS": "pos"}), on='pos', how='inner')
            clear_dfs(df_fst_combined, df_fst_merged, df_fst_final)
            print(f"FST files for {pop} chromosome {chr_num} merged successfully")
            print(df_fst_final.head(2))
            preview_output(df_empirical, "adding mean_fst")

        # Process NSL files
        nsl_file = os.path.join(normed_path, f"norm_nsl_{pop}_{chr_num}.tsv")
        if os.path.exists(nsl_file):
            print(f"Reading NSL file for {pop} chromosome {chr_num}")
            df_nsl = pd.read_csv(nsl_file, sep='\t').drop_duplicates(subset='pos').rename(columns={'rs': 'rsid_nsl'}).astype({'pos': 'Int64'})
            df_empirical = pd.merge(df_empirical, df_nsl[['pos', 'rsid_nsl', 'norm_nsl']], on='pos', how='outer')
            clear_dfs(df_nsl)
            print(f"NSL file for {pop} chromosome {chr_num} read and merged successfully")
            print(f"Word count for NSL file: {len(df_nsl)}")
            print(df_nsl.head(2))
            preview_output(df_empirical, "adding norm_nsl")

        # Process IHH12 files
        ihh12_file = os.path.join(normed_path, f"norm_ihh12_{pop}_{chr_num}.tsv")
        if os.path.exists(ihh12_file):
            print(f"Reading IHH12 file for {pop} chromosome {chr_num}")
            df_ihh12 = pd.read_csv(ihh12_file, sep='\t').drop_duplicates(subset='pos').rename(columns={'id': 'rsid_ihh12'}).astype({'pos': 'Int64'})
            df_empirical = pd.merge(df_empirical, df_ihh12[['pos', 'rsid_ihh12', 'norm_ihh12']], on='pos', how='outer')
            clear_dfs(df_ihh12)
            print(f"IHH12 file for {pop} chromosome {chr_num} read and merged successfully")
            print(f"Word count for IHH12 file: {len(df_ihh12)}")
            print(df_ihh12.head(2))
            preview_output(df_empirical, "adding norm_ihh12")

        # Create inclusive rsid column
        df_empirical['rsid_inclusive'] = df_empirical[['rsid_nsl', 'rsid_ihh12']].bfill(axis=1).iloc[:, 0]

        # Remove duplicates by rsid_inclusive and pos
        df_empirical = df_empirical.drop_duplicates(subset=['rsid_inclusive', 'pos'])
        preview_output(df_empirical, "removing duplicates by rsid_inclusive and pos")

        # Process DELIHH files
        delihh_file = os.path.join(normed_path, f"norm_delihh_{pop}_{chr_num}.tsv")
        if os.path.exists(delihh_file):
            print(f"Reading DELIHH file for {pop} chromosome {chr_num}")
            df_delihh = pd.read_csv(delihh_file, sep='\t').drop_duplicates(subset='Index').rename(columns={'ID': 'rsid_delihh'}).astype({'Index': 'Int64'})
            df_empirical = pd.merge(df_empirical, df_delihh[['rsid_delihh', 'norm_iHS', 'norm_delihh']], left_on='rsid_inclusive', right_on='rsid_delihh', how='outer')
            df_empirical.drop(columns=['rsid_delihh'], inplace=True, errors='ignore')
            clear_dfs(df_delihh)
            print(f"DELIHH file for {pop} chromosome {chr_num} read and merged successfully")
            print(f"Word count for DELIHH file: {len(df_delihh)}")
            print(df_delihh.head(2))
            preview_output(df_empirical, "adding norm_delihh")
        
        df_empirical.to_csv(empirical_file, sep='\t', index=False)
        print(f"Empirical file for {pop} chromosome {chr_num} saved successfully")
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
                    print(f"Reading XPEHH file {i} for {pop} chromosome {chr_num}")
                    df_xpehh = pd.read_csv(xpehh_file, sep='\t')
                    df_xpehh = df_xpehh.drop_duplicates(subset='Index').rename(columns={'ID': 'rsid', 'norm_XPEHH': f'norm_XPEHH_{i}'}).astype({'Index': 'Int64'})
                    df_xpehh_combined.append(df_xpehh)
                    print(f"XPEHH file {i} for {pop} chromosome {chr_num} read successfully")
                    print(f"Word count for XPEHH file {i}: {len(df_xpehh)}")
                    print(df_xpehh.head(2))

            # Use inner join for calculating max_xpehh
            print(f"Merging XPEHH files for {pop} chromosome {chr_num}")
            df_xpehh_merged = reduce(lambda left, right: pd.merge(left, right, on='rsid', how='inner'), df_xpehh_combined)
            print(f"Temporary merged XPEHH file for {pop} chromosome {chr_num}")
            print(df_xpehh_merged.head(2))
            xpehh_columns = [col for col in df_xpehh_merged.columns if col.startswith('norm_XPEHH_')]
            df_xpehh_merged['max_xpehh'] = df_xpehh_merged[xpehh_columns].apply(lambda row: row.max(), axis=1)

            df_empirical = pd.merge(df_empirical, df_xpehh_merged[['Index', 'rsid', 'max_xpehh']], left_on='rsid_inclusive', right_on='rsid', how='outer')
            df_empirical.drop(columns=['Index', 'rsid', 'rsid_nsl', 'rsid_ihh12'], inplace=True, errors='ignore')
            df_empirical.rename(columns={'rsid_inclusive': 'rsid'}, inplace=True)
            df_empirical.to_csv(empirical_file, sep='\t', index=False)
            clear_dfs(df_xpehh_combined, df_xpehh_merged)
            print(f"XPEHH files for {pop} chromosome {chr_num} merged successfully")
            preview_output(df_empirical, "adding max_xpehh")

        # Convert mean_fst to 5 significant figures, replace NaNs with pd.NA, reorder rows and columns
        df_empirical['mean_fst'] = df_empirical['mean_fst'].apply(lambda x: round(x, 5))
        df_empirical = df_empirical.sort_values(by='pos')
        df_empirical = df_empirical.fillna(pd.NA)
        df_empirical = df_empirical[['rsid'] + [col for col in df_empirical.columns if col != 'rsid']]

        # Remove duplicate pos rows
        df_empirical = df_empirical.drop_duplicates(subset='pos')
        preview_output(df_empirical, "removing duplicate pos rows")

        # Save final DataFrame
        df_empirical.to_csv(empirical_file, sep='\t', index=False)
        clear_dfs(df_empirical)
        gc.collect()
        print(f"Final empirical file for {pop} chromosome {chr_num} saved successfully")

# Process each chromosome sequentially
for chr_num in range(1, 23):
    process_chromosome(chr_num)

print("Processing completed.")
