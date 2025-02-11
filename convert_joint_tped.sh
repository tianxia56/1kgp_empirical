#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000

# List of populations to run vcftools
target_pops=("YRI" "CEU" "CHB" "BEB")

# Function to create joint map file
create_joint_map() {
  local chr=$1
  local joint_map="/home/tx56/ycga_work/1kgp_joint/joint.map.chr${chr}.map"

  # Combine map files and keep only rows that appear 4 times
  cat /home/tx56/ycga_work/1kgp_map/*chr${chr}.map | awk '{print $2}' | sort | uniq -c | awk '$1 == 4 {print $2}' > $joint_map
}

# Run plink for each population and chromosome
for pop in ${target_pops[@]}
do
  for i in {1..22}
  do
    # Create joint map file
    create_joint_map $i

    # Convert plink file to tped
    plink --file /home/tx56/ycga_work/1kgp_map/${pop}.chr${i} --recode 01 transpose --cm-map /home/tx56/selscan_maps/${pop}/${pop}_map_chr_${i}.txt ${i} --biallelic-only --exclude /home/tx56/1KGP/exclude_rsid.txt --output-missing-genotype 1 --extract /home/tx56/ycga_work/1kgp_joint/joint.map.chr${i}.map --out /home/tx56/ycga_work/1kgp_joint/${pop}.${i}.joint

  done
done
