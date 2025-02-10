#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000

# List of populations to run selscan
target_pops=("YRI" "CEU" "CHB" "BEB")

# Function to compare and switch alleles if necessary
compare_and_switch_alleles() {
  local pop=$1
  local chr=$2
  local tped_bim_file="/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped/${pop}.${chr}.AA.bim"
  local bim_file="/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_bfile/${pop}.${chr}.AA.bim"
  local temp_bim_file="/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped/temp_${pop}.${chr}.AA.bim"

  # Compare .bim files and switch alleles if necessary
  awk 'NR==FNR{a[$2]=$2" "$5" "$6; next}
      {if($5!=a[$2] || $6!=a[$3]) print $1,$2,$3,$4,a[$2],a[$3]; else print $0}' $bim_file $tped_bim_file > $temp_bim_file

  # Replace the original .bim file with the modified one
  mv $temp_bim_file $tped_bim_file
}

mkdir -p /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped
# Loop through each population and chromosome
for pop in ${target_pops[@]}
do
  for i in {1..22}
  do
    # Compare and switch alleles if necessary
    compare_and_switch_alleles $pop $i

    # Convert plink file to tped
    plink --bfile /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped/${pop}.${i}.AA --recode 01 transpose --cm-map /home/tx56/selscan_maps/${pop}/${pop}_map_chr_${i}.txt ${i} --biallelic-only --exclude /home/tx56/1KGP/exclude_rsid.txt --output-missing-genotype 1 --out /home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped/${pop}.${i}.AA
  done
done

