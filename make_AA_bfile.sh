#!/bin/bash

mkdir -p /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_tped

# List of populations to run selscan
target_pops=("YRI" "CEU" "CHB" "BEB")

# Loop through each population and chromosome
for pop in ${target_pops[@]}
do
  for i in {1..22}
  do
    plink --file /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_map/${pop}.chr${i} --make-bed --extract /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_bfile/${pop}.${i}.AA.bim --out /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_tped/${pop}.${i}.AA
  done
done


