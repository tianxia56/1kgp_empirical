import os

def extract_and_clean_columns(input_file, output_file_map, output_file_hap):
    with open(input_file, 'r') as file, open(output_file_map, 'w') as map_file, open(output_file_hap, 'w') as hap_file:
        for line in file:
            columns = line.split()
            if len(columns) < 5:
                continue
            map_file.write(' '.join(columns[:4]) + '\n')  # Extract the first four columns with whitespace separator
            cleaned_columns = [col.strip() for col in columns[4:]]  # Clean columns
            hap_file.write(' '.join(cleaned_columns) + '\n')  # Write cleaned columns with whitespace separator

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

def main():
    target_pops = ["YRI", "CEU", "CHB", "BEB"]
    path = "/home/tx56/ycga_work/1kgp_tped"

    for pop in target_pops:
        for i in range(1, 23):
            input_file = f"{path}/{pop}.{i}.tped"
            
            if not os.path.exists(input_file):
                print(f"Skipping {pop} chromosome {i} as {input_file} does not exist.")
                continue

            os.makedirs("hapbin", exist_ok=True)
            output_file_hap = f"/home/tx56/ycga_work/1kgp_ihs/{pop}.{i}.hap"
            output_file_map = f"/home/tx56/ycga_work/1kgp_ihs/{pop}.{i}.map"
            extract_and_clean_columns(input_file, output_file_map, output_file_hap)
            print(f"Extracted columns saved to {output_file_map} and {output_file_hap}")
            create_map_file(output_file_map)

if __name__ == "__main__":
    main()
