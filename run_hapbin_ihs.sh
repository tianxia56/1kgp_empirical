#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000

# Run the Python script to convert TPED to HAP
python ihs_tped_to_hap.py

# List of populations to run selscan
target_pops=("YRI" "CEU" "CHB" "BEB")

# Loop through each population and chromosome
for pop in ${target_pops[@]}
do
  for i in {1..22}
  do
    /home/tx56/hapbin/build/ihsbin --hap /home/tx56/ycga_work/1kgp_ihs/${pop}.${i}.hap --map /home/tx56/ycga_work/1kgp_ihs/${pop}.${i}.map --out /home/tx56/ycga_work/1kgp_ihs/${pop}.${i}.ihs
  done
done
