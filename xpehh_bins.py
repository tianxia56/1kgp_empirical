import os
import pandas as pd
import numpy as np
from collections import defaultdict

# Function to read and process xpehh files
def process_xpehh(file_path):
    df = pd.read_csv(file_path, sep='\t', header=0, names=['Index', 'ID', 'Freq', 'iHH_A1', 'iHH_B1', 'iHH_P1', 'XPEHH', 'std_XPEHH'])
    df['Freq'] = pd.to_numeric(df['Freq'], errors='coerce')  # Ensure 'Freq' is float
    return df[['Index', 'Freq', 'XPEHH']].rename(columns={'Freq': 'daf'})

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
xpehh_dir = os.path.join(base_path, '1kgp_xpehh')

# Get list of all xpehh files and group by population pair
xpehh_files = [f for f in os.listdir(xpehh_dir) if f.endswith('.joint.xpehh')]
pop_pair_files = defaultdict(list)

for file in xpehh_files:
    pop_pair = file.split('.')[0]  # Extract population pair from file name
    pop_pair_files[pop_pair].append(file)

# Process each population pair
for pop_pair, files in pop_pair_files.items():
    xpehh_data = pd.DataFrame()

    # Process all chromosome files for each population pair
    for file in files:
        file_path = os.path.join(xpehh_dir, file)
        if os.path.exists(file_path):
            xpehh_data = pd.concat([xpehh_data, process_xpehh(file_path)], ignore_index=True)

    # Debugging: Check if the dataframe has the expected columns
    print(f"Processing population pair {pop_pair}")
    print(f"xpehh_data columns: {xpehh_data.columns}")

    # Drop rows with NaN values in the daf column
    if 'daf' in xpehh_data.columns:
        xpehh_data.dropna(subset=['daf'], inplace=True)

    # Bin data and calculate mean and std based on the distribution of daf
    output_dir = os.path.join(base_path, 'bins')
    os.makedirs(output_dir, exist_ok=True)
    
    if not xpehh_data.empty:
        xpehh_binned = bin_data(xpehh_data, 'XPEHH')
        xpehh_binned.to_csv(f'{output_dir}/{pop_pair}_xpehh_bin.csv', index=False)

    print(f"Binned data saved to the '{output_dir}' directory for population pair {pop_pair}.")
