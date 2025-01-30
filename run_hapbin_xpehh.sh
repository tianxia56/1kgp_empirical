#!/bin/bash
#SBATCH --partition=week
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50000

# Run the Python script to convert TPED to HAP
python xpehh_tped_to_hap.py

# List of populations to run selscan
target_pops=("YRI" "CEU" "CHB" "BEB")

# Loop through each pair of populations and chromosome
for pop_a in ${target_pops[@]}
do
  for pop_b in ${target_pops[@]}
  do
    if [ "$pop_a" != "$pop_b" ]; then
      for i in {1..22}
      do
        /home/tx56/hapbin/build/xpehhbin --hapA /home/tx56/ycga_work/1kgp_joint/${pop_a}.${i}.joint.hap --hapB /home/tx56/ycga_work/1kgp_joint/${pop_b}.${i}.joint.hap --map /home/tx56/ycga_work/1kgp_joint/${pop_a}.${i}.joint.map --out /home/tx56/ycga_work/1kgp_xpehh/${pop_a}_${pop_b}.${i}.joint.xpehh
      done
    fi
  done
done
