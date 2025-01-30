#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000

# List of populations to run selscan
target_pops=("YRI" "CEU" "CHB" "BEB")

# Loop through each population and chromosome
for pop in ${target_pops[@]}
do
  for i in {1..22}
  do
    plink --file /home/tx56/ycga_work/1kgp_map/${pop}.chr${i} --recode 01 transpose --cm-map /home/tx56/selscan_maps/${pop}/${pop}_map_chr_${i}.txt ${i} --biallelic-only --exclude /home/tx56/1KGP/exclude_rsid.txt --output-missing-genotype 1 --out /home/tx56/ycga_work/1kgp_tped/${pop}.${i}
  done
done
