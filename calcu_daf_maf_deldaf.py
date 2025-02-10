import pandas as pd
import os

# Define target populations and chromosomes
target_pops = ["YRI", "CEU", "CHB", "BEB"]
chromosomes = range(1, 23)

# Function to safely read and inspect the file
def read_tped_file(file_path):
    try:
        df = pd.read_csv(file_path, sep=' ', header=None)
        print(f"Read {file_path} successful with shape {df.shape}")
        return df
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return pd.DataFrame()

# Loop through each chromosome
for chromosome in chromosomes:
    # Initialize list to store data for merging
    data_list = []

    # Loop through each population to compute DAFs
    for pop in target_pops:
        # Define file path
        tped_file = f"/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped/{pop}.{chromosome}.AA.tped"
        
        # Read TPED file
        pop_data = read_tped_file(tped_file)

        if pop_data.empty:
            print(f"Skipping {tped_file} due to read error or empty file.")
            continue

        # Verify the file contains enough columns
        if pop_data.shape[1] < 5:
            print(f"{tped_file} does not contain enough columns.")
            continue

        # Extract positions and alleles
        pop_data = pop_data.drop_duplicates(subset=[3])
        positions = pop_data.iloc[:, 3].astype(int)
        alleles = pop_data.iloc[:, 4:]

        # Compute derived allele frequency (DAF)
        daf = alleles.apply(lambda x: (x == 1).mean(), axis=1).round(4)

        # Store the data in the list
        daf_df = pd.DataFrame({'pos': positions, f'daf_{pop}': daf})
        data_list.append(daf_df)

        # Debug: Check for duplicates in positions
        if daf_df['pos'].duplicated().any():
            print(f"Warning: Duplicates found in positions for {pop} on chromosome {chromosome}")
            print(daf_df['pos'][daf_df['pos'].duplicated()])

    if not data_list:
        print(f"No data available for chromosome {chromosome}. Skipping...")
        continue

    # Merge all populations' data by position to compute deldaf
    merged_data = data_list[0]
    for df in data_list[1:]:
        merged_data = merged_data.merge(df, on='pos', how='inner')

    # Compute ΔDAF for each SNP where no data is missing
    daf_columns = [col for col in merged_data.columns if col.startswith('daf_')]
    merged_data['deldaf'] = merged_data[daf_columns].apply(lambda x: round(x.max() - x.min(), 4) if x.notna().all() else None, axis=1)

    # Prepare output for each population
    for pop in target_pops:
        # Define file paths
        output_file = f"/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_daf/{pop}.{chromosome}.AA.tsv"

        # Extract the specific DAF column for the current population
        daf_column = f'daf_{pop}'
        if daf_column not in merged_data.columns:
            print(f"DAF column {daf_column} not found in merged data for {pop} on chromosome {chromosome}")
            continue

        # Prepare the results table by performing an outer join with the `deldaf` table
        pop_daf_df = data_list[target_pops.index(pop)]
        results = pop_daf_df.merge(merged_data[['pos', 'deldaf']], on='pos', how='outer')
        results.rename(columns={daf_column: 'daf'}, inplace=True)

        # Re-read TPED file to compute MAF
        pop_data = read_tped_file(f"/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped/{pop}.{chromosome}.AA.tped")

        if pop_data.empty:
            print(f"Skipping {tped_file} for final output due to read error or empty file.")
            continue

        # Extract positions and alleles
        pop_data = pop_data.drop_duplicates(subset=[3])
        positions = pop_data.iloc[:, 3].astype(int)
        alleles = pop_data.iloc[:, 4:]

        # Compute derived allele frequency (DAF)
        daf = alleles.apply(lambda x: (x == 1).mean(), axis=1).round(4)

        # Compute minor allele frequency (MAF)
        maf = daf.apply(lambda x: round(min(x, 1 - x), 4))
        
        # Merge MAF with the results
        maf_df = pd.DataFrame({'pos': positions, 'maf': maf})
        results = results.merge(maf_df, on='pos', how='left')
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Write results to a file, ensuring pos is written as an integer
        results.to_csv(output_file, sep='\t', index=False)

