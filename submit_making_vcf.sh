#!/bin/bash

# List of populations to extract
target_pops=("YRI" "CEU" "CHB" "BEB")

# Create SLURM job scripts for each population
for pop in ${target_pops[@]}
do
  job_script="/home/tx56/ycga_work/1kgp_vcf/extract_${pop}.sh"
  cat <<EOT > $job_script
#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000

# Extract VCF files for population $pop
for i in {1..22}
do
  vcftools --gzvcf /home/tx56/1KGP/ALL.chr\${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
           --keep /home/tx56/ycga_work/empirical/pop_files/${pop}.txt \
           --recode \
           --recode-INFO-all \
           --maf 0.05 \
           --min-alleles 2 \
           --max-alleles 2 \
           --out /home/tx56/ycga_work/1kgp_vcf/${pop}.chr\${i}
done
EOT
  # Submit the job
  sbatch $job_script
done
