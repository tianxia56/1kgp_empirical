# Load necessary library
library(data.table)

# Define target populations and chromosomes
target_pops <- c("YRI", "CEU", "CHB", "BEB")
chromosomes <- 1:22

# Loop through each chromosome
for (i in chromosomes) {
  # Initialize a list to store DAFs for each population
  daf_list <- list()
  pos_list <- list()
  
  # Loop through each population to compute DAFs
  for (pop in target_pops) {
    # Define file paths
    tped_file <- sprintf("/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped/%s.%d.AA.tped", pop, i)
    
    # Read TPED file
    pop_data <- fread(tped_file, select = 5:ncol(fread(tped_file)))
    
    # Compute derived allele frequency (DAF)
    daf <- rowMeans(pop_data == 1, na.rm = TRUE)
    
    # Store DAFs and positions in the list
    daf_list[[pop]] <- daf
    pos_list[[pop]] <- fread(tped_file, select = 4)[[1]]
  }
  
  # Combine DAFs and positions into a data table
  daf_combined <- do.call(cbind, daf_list)
  pos_combined <- do.call(cbind, pos_list)
  
  # Compute ΔDAF for each SNP
  deldaf <- apply(daf_combined, 1, function(x) {
    if (sum(!is.na(x)) < length(target_pops)) {
      return(NA)
    } else {
      return(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    }
  })
  
  # Loop through each population to compute MAF and save results
  for (pop in target_pops) {
    # Define file paths
    tped_file <- sprintf("/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_AA_tped/%s.%d.AA.tped", pop, i)
    output_file <- sprintf("/home/tx56/palmer_scratch/deepsweep_empirical/1kgp_daf/%s.%d.AA.tsv", pop, i)
    
    # Create output directory if it doesn't exist
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    
    # Read TPED file
    pop_data <- fread(tped_file, select = 5:ncol(fread(tped_file)))
    
    # Compute derived allele frequency (DAF)
    daf <- rowMeans(pop_data == 1, na.rm = TRUE)
    
    # Compute minor allele frequency (MAF)
    maf <- round(pmin(daf, 1 - daf), 4)
    
    # Create results table
    results <- data.table(pos = pos_combined[, pop], daf = round(daf, 4), deldaf = round(deldaf, 4), maf = maf)
    
    # Write results to a file
    fwrite(results, output_file, sep = "\t", quote = FALSE)
  }
}
