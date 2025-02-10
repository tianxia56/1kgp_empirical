import os
import pandas as pd
import numpy as np

# Function to read and process nsl files
def process_nsl(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None, names=['ID', 'pos', 'p1', 'sl1', 'sl0', 'nsl'], skiprows=1)
    df['p1'] = pd.to_numeric(df['p1'], errors='coerce')  # Ensure 'p1' is float
    return df[['pos', 'p1', 'nsl']].rename(columns={'p1': 'daf'})

# Function to read and process ihs files
def process_ihs(file_path):
    df = pd.read_csv(file_path, sep='\t', header=0, names=['Index', 'ID', 'Freq', 'iHH_0', 'iHH_1', 'iHS', 'Std_iHS'])
    df['Freq'] = pd.to_numeric(df['Freq'], errors='coerce')  # Ensure 'Freq' is float
    df['delihh'] = df['iHH_1'] - df['iHH_0']  # Compute delihh
    return df[['Index', 'Freq', 'iHS']].rename(columns={'Freq': 'daf'}), df[['Index', 'Freq', 'delihh']].rename(columns={'Freq': 'daf'})

# Function to read and process ihh12 files
def process_ihh12(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None, names=['id', 'pos', 'p1', 'ihh12'], skiprows=1)
    df['p1'] = pd.to_numeric(df['p1'], errors='coerce')  # Ensure 'p1' is float
    return df[['pos', 'p1', 'ihh12']].rename(columns={'p1': 'daf'})

# Function to bin data and calculate mean and std based on the distribution of daf
def bin_data(df, value_col, num_bins=20):
    # Sort the dataframe by daf
    df = df.sort_values(by='daf')
    
    # Calculate bin edges based on the distribution of daf
    bin_edges = np.linspace(df['daf'].min(), df['daf'].max(), num_bins + 1)
    
    # Assign each row to a bin
    df['bin'] = pd.cut(df['daf'], bins=bin_edges, include_lowest=True)
    
    # Calculate mean and std for each bin
    binned = df.groupby('bin', observed=True)[value_col].agg(['mean', 'std']).reset_index()
    
    # Extract the right edge of each bin for labeling
    binned['bin'] = binned['bin'].apply(lambda x: x.right)
    
    return binned

# Paths to the directories containing the files
base_path = '/home/tx56/palmer_scratch/deepsweep_empirical'
nsl_dir = os.path.join(base_path, '1kgp_nsl')
ihs_dir = os.path.join(base_path, '1kgp_ihs')
ihh12_dir = os.path.join(base_path, '1kgp_ihh12')

# List of populations
populations = ['BEB', 'CEU', 'CHB', 'JPT', 'YRI']  # Add other populations as needed

# Process each population separately
for pop in populations:
    # Initialize empty dataframes for concatenation
    nsl_data = pd.DataFrame()
    ihs_data = pd.DataFrame()
    delihh_data = pd.DataFrame()
    ihh12_data = pd.DataFrame()

    # Process all chromosome files for each statistic
    for chr_num in range(1, 23):  # Assuming chromosomes 1 to 22
        nsl_file = os.path.join(nsl_dir, f'{pop}.{chr_num}.AA.nsl.out')
        ihs_file = os.path.join(ihs_dir, f'{pop}.{chr_num}.ihs')
        ihh12_file = os.path.join(ihh12_dir, f'{pop}.{chr_num}.AA.ihh12.out')

        if os.path.exists(nsl_file):
            nsl_data = pd.concat([nsl_data, process_nsl(nsl_file)], ignore_index=True)
        
        if os.path.exists(ihs_file):
            ihs_df, delihh_df = process_ihs(ihs_file)
            ihs_data = pd.concat([ihs_data, ihs_df], ignore_index=True)
            delihh_data = pd.concat([delihh_data, delihh_df], ignore_index=True)
        
        if os.path.exists(ihh12_file):
            ihh12_data = pd.concat([ihh12_data, process_ihh12(ihh12_file)], ignore_index=True)

    # Debugging: Check if the dataframes have the expected columns
    print(f"Processing population {pop}")
    print(f"nsl_data columns: {nsl_data.columns}")
    print(f"ihs_data columns: {ihs_data.columns}")
    print(f"delihh_data columns: {delihh_data.columns}")
    print(f"ihh12_data columns: {ihh12_data.columns}")

    # Drop rows with NaN values in the daf column
    if 'daf' in nsl_data.columns:
        nsl_data.dropna(subset=['daf'], inplace=True)
    if 'daf' in ihs_data.columns:
        ihs_data.dropna(subset=['daf'], inplace=True)
    if 'daf' in delihh_data.columns:
        delihh_data.dropna(subset=['daf'], inplace=True)
    if 'daf' in ihh12_data.columns:
        ihh12_data.dropna(subset=['daf'], inplace=True)

    # Bin data and calculate mean and std based on the distribution of daf for each dataset separately
    output_dir = os.path.join(base_path, 'bins')
    os.makedirs(output_dir, exist_ok=True)
    
    if not nsl_data.empty:
        nsl_binned = bin_data(nsl_data, 'nsl')
        nsl_binned.to_csv(f'{output_dir}/{pop}_nsl_bin.csv', index=False)
    if not ihs_data.empty:
        ihs_binned = bin_data(ihs_data, 'iHS')
        ihs_binned.to_csv(f'{output_dir}/{pop}_ihs_bin.csv', index=False)
    if not delihh_data.empty:
        delihh_binned = bin_data(delihh_data, 'delihh')
        delihh_binned.to_csv(f'{output_dir}/{pop}_delihh_bin.csv', index=False)
    if not ihh12_data.empty:
        ihh12_binned = bin_data(ihh12_data, 'ihh12')
        ihh12_binned.to_csv(f'{output_dir}/{pop}_ihh12_bin.csv', index=False)

    print(f"Binned data saved to the '{output_dir}' directory for population {pop}.")
