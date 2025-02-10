#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100000

#bash submit_making_vcf.sh
#bash submit_plink_map.sh
#bash submit_fst.sh
#bash submit_making_bfile.sh

#Rscript convert_AA_bfile.R
#bash AA_bfile.sh
#bash convert_ihs_AA_tped.sh
#bash convert_joint_AA.tped.sh

#python ihs_tped_to_hap.py
#bash run_hapbin_ihs.sh 
#python xpehh_tped_to_hap.py
#bash run_hapbin_xpehh.sh
#bash run_selscan_AA_nsl_ihh12.sh

#python calcu_daf_maf_deldaf.py
#python one_pop_bins.py
#python xpehh_bins.py
#Rscript norm.R

#python map_index_outputs.py

#make_isafe_haps
#add_isafe.py
#plot_components.py

#plot_chr.py
