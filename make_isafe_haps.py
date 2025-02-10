import os
import re

def extract_and_clean_columns(input_file, output_file_map, output_file_hap, sel_region):
    with open(input_file, 'r') as file, open(output_file_map, 'w') as map_file, open(output_file_hap, 'w') as hap_file:
        for line in file:
            columns = line.split()
            if len(columns) < 5:
                continue
            pos = int(columns[3])  # 4th column is the position
            if sel_region[0] <= pos <= sel_region[1]:
                map_file.write(' '.join(columns[:4]) + '\n')  # Extract the first four columns with whitespace separator
                cleaned_columns = [col.strip() for col in columns[4:]]  # Clean columns
                hap_file.write(f"{pos}\t" + '\t'.join(cleaned_columns) + '\n')  # Write position as first column and cleaned columns with tab separator

def create_map_file(original_map_file):
    temp_file = original_map_file + ".tmp"
    with open(original_map_file, 'r') as map_file, open(temp_file, 'w') as new_file:
        for line in map_file:
            columns = line.split()
            # Ensure the columns are correctly formatted for ASCII and convert scientific notation to numeric, except for the second column
            formatted_columns = [col.encode('ascii', 'ignore').decode('ascii') for col in columns]
            formatted_columns = [
                f"{int(col)}" if idx in [0, 3] else col if idx == 1 else f"{float(col):.6f}" if col.replace('.', '', 1).isdigit() else col
                for idx, col in enumerate(formatted_columns)
            ]
            # Join columns with a single space and remove any non-printable characters
            cleaned_line = ' '.join(formatted_columns).strip()
            new_file.write(cleaned_line + '\n')
    os.replace(temp_file, original_map_file)
    print(f"Map file created: {original_map_file}")

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

def process_input_definition(input_definition):
    gene, chr, sel_region, pop = parse_input_definition(input_definition)
    path = "/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped"
    input_file = f"{path}/{pop}.{chr}.AA.tped"
    
    if not os.path.exists(input_file):
        print(f"Skipping {pop} chromosome {chr} as {input_file} does not exist.")
        return

    os.makedirs("/home/tx56/palmer_scratch/deepsweep_empirical/sel_loci", exist_ok=True)
    output_file_hap = f"/home/tx56/palmer_scratch/deepsweep_empirical/sel_loci/{gene}.{pop}.{chr}.hap"
    output_file_map = f"/home/tx56/palmer_scratch/deepsweep_empirical/sel_loci/{gene}.{pop}.{chr}.map"
    extract_and_clean_columns(input_file, output_file_map, output_file_hap, sel_region)
    print(f"Extracted columns saved to {output_file_map} and {output_file_hap}")
    create_map_file(output_file_map)

    # Ensure positions in hap file are integers and sorted
    with open(output_file_hap, 'r') as hap_file:
        lines = hap_file.readlines()
    lines.sort(key=lambda x: int(x.split()[0]))  # Sort lines by the first column (position)
    with open(output_file_hap, 'w') as hap_file:
        for line in lines:
            columns = line.split()
            columns[0] = str(int(columns[0]))  # Ensure the first column (position) is an integer
            hap_file.write('\t'.join(columns) + '\n')

    # Run isafe command
    os.system(f"isafe --input {output_file_hap} --output ../sel_loci/{gene}.{pop} --format hap --IgnoreGaps")

def main():
    input_definitions = [
        "LCT,2:134070080..139070080,CEU",
        "SLC24A5,15:45924019..50924019,BEB",
        "EDAR,2:106411000..111411000,CHB",
        "DARC,1:156677251..161677251,YRI",
        "SLC45A2,5:31451691..36451691,CEU"  # Centered SLC45A2 gene and variant with a 5MB window
    ]
    
    for input_definition in input_definitions:
        process_input_definition(input_definition)

if __name__ == "__main__":
    main()
