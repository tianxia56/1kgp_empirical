#!/bin/bash

# List of populations to run selscan
target_pops=("YRI" "CEU" "CHB" "BEB")

# Create SLURM job scripts for each population
for pop in ${target_pops[@]}
do
  job_script="/home/tx56/ycga_work/1kgp_ihh12/run_selscan_${pop}_ihh12.sh"
  cat <<EOT > $job_script
#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8000

# Run selscan for population $pop
for i in {1..22}
do
  selscan --ihh12 \
          --vcf /home/tx56/ycga_work/1kgp_vcf/${pop}.chr\${i}.recode.vcf \
          --map /home/tx56/ycga_work/1kgp_map/${pop}.chr\${i}.map \
          --out /home/tx56/ycga_work/1kgp_ihh12/${pop}.chr\${i} \
          --threads 8
done
EOT
  # Submit the job
  sbatch $job_script
done
