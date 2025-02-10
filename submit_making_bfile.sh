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
    plink --file /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_map/${pop}.chr${i} --make-bed --biallelic-only --exclude /home/tx56/1KGP/exclude_rsid.txt --out /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_bfile/${pop}.${i}
  done
done
