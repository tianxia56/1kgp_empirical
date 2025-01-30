#!/bin/bash

# List of populations to convert
target_pops=("YRI" "CEU" "CHB" "BEB")

# Create SLURM job scripts for each population
for pop in ${target_pops[@]}
do
  job_script="/home/tx56/ycga_work/1kgp_map/convert_${pop}.sh"
  cat <<EOT > $job_script
#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000

# Convert VCF files to PLINK map files for population $pop
for i in {1..22}
do
  plink --vcf /home/tx56/ycga_work/1kgp_vcf/${pop}.chr\${i}.recode.vcf \
        --recode \
        --out /home/tx56/ycga_work/1kgp_map/${pop}.chr\${i}
done
EOT
  # Submit the job
  sbatch $job_script
done
