#!/bin/bash
#SBATCH --partition=week
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8000

target_pops=("YRI" "CEU" "CHB" "BEB")
for pop in ${target_pops[@]}
do
  for i in {1..22}
  do
    tped_file="/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped/${pop}.${i}.AA.tped"
    
    mkdir -p /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_nsl
    mkdir -p /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_ihh12
    
    selscan --nsl --tped $tped_file --out /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_nsl/${pop}.${i}.AA --threads 8
    selscan --ihh12 --tped $tped_file --out /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_ihh12/${pop}.${i}.AA --threads 8
  done
done
