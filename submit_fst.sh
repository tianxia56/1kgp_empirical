#!/bin/bash

# Define the panel file
panel_file="/home/tx56/1KGP/integrated_call_samples_v3.20130502.ALL.panel"

# Create a directory to store the population files
mkdir -p /home/tx56/ycga_work/empirical/pop_files

# Extract unique populations
populations=$(cut -f2 $panel_file | tail -n +2 | sort | uniq)

# Generate .txt files for each population
for pop in $populations
do
  grep -w $pop $panel_file | cut -f1 > /home/tx56/ycga_work/empirical/pop_files/${pop}.txt
done

# List of populations to compare
target_pops=("YRI" "CEU" "CHB" "BEB")

# Create SLURM job scripts for each pairwise comparison
job_count=0
for i in ${!target_pops[@]}
do
  for j in $(seq $((i + 1)) ${#target_pops[@]})
  do
    if [ $j -lt ${#target_pops[@]} ]; then
      pop_a=${target_pops[$i]}
      pop_b=${target_pops[$j]}
      job_count=$((job_count + 1))
      job_script="/home/tx56/ycga_work/empirical/job_${job_count}.sh"
      cat <<EOT > $job_script
#!/bin/bash
#SBATCH --partition=week
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100000

# Run VCFtools for each VCF file
for vcf_file in /home/tx56/1KGP/ALL.chr{1..22}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
do
  vcftools --gzvcf \$vcf_file \
           --weir-fst-pop /home/tx56/ycga_work/empirical/pop_files/${pop_a}.txt \
           --weir-fst-pop /home/tx56/ycga_work/empirical/pop_files/${pop_b}.txt \
           --out /home/tx56/ycga_work/1kgp_fst/${pop_a}.${pop_b}.\$(basename \$vcf_file .vcf.gz)
done
EOT
      # Submit the job
      sbatch $job_script
    fi
  done
done
