import os
import re
import pandas as pd

def parse_input_definition(input_definition):
    match = re.match(r"(\w+),(\d+):(\d+)\.\.(\d+),(\w+)", input_definition)
    if match:
        gene = match.group(1)
        chr = int(match.group(2))
        start = int(match.group(3))
        end = int(match.group(4))
        pop = match.group(5)
        return gene, chr, (start, end), pop
    else:
        raise ValueError("Input definition format is incorrect")

def extract_empirical_data(gene, chr, sel_region, pop):
    empirical_file = f"../outputs/{pop}.{chr}.empirical.tsv"
    if not os.path.exists(empirical_file):
        print(f"Empirical file {empirical_file} does not exist.")
        return None

    empirical_data = pd.read_csv(empirical_file, sep='\t')
    filtered_data = empirical_data[(empirical_data['chr'] == chr) & 
                                   (empirical_data['pos'] >= sel_region[0]) & 
                                   (empirical_data['pos'] <= sel_region[1])]
    return filtered_data

def merge_with_isafe(gene, pop, empirical_data):
    isafe_file = f"../sel_loci/{gene}.{pop}.iSAFE.out"
    if not os.path.exists(isafe_file):
        print(f"iSAFE file {isafe_file} does not exist.")
        return empirical_data

    isafe_data = pd.read_csv(isafe_file, sep='\t')
    merged_data = pd.merge(empirical_data, isafe_data[['POS', 'iSAFE']], how='left', left_on='pos', right_on='POS')
    merged_data.drop(columns=['POS'], inplace=True)  # Drop the POS column after merging
    return merged_data

def process_meta_info(input_definitions):
    for input_definition in input_definitions:
        gene, chr, sel_region, pop = parse_input_definition(input_definition)

        # Extract empirical data
        empirical_data = extract_empirical_data(gene, chr, sel_region, pop)
        if empirical_data is not None:
            # Merge with iSAFE data
            merged_data = merge_with_isafe(gene, pop, empirical_data)
            output_stats_file = f"../sel_loci/{gene}.{pop}.stats"
            merged_data.to_csv(output_stats_file, sep='\t', index=False)
            print(f"Stats saved to {output_stats_file}")

if __name__ == "__main__":
    input_definitions = [
        "LCT,2:134070080..139070080,CEU",
        "SLC24A5,15:45924019..50924019,BEB",  # Changed to BEB
        "EDAR,2:106411000..111411000,CHB",
        "DARC,1:156677251..161677251,YRI",
        "SLC45A2,5:31451691..36451691,CEU"  # Centered SLC45A2 variant with a 5MB window
    ]
    process_meta_info(input_definitions)
