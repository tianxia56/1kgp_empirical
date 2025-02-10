# Create the output directory if it doesn't exist
output_dir <- "/home/tx56/palmer_scratch/deepsweep_empirical/normed"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to read and normalize data
normalize_data <- function(input_file, norm_file, output_file, value_col, col_names, required_columns) {
  cat("Processing file:", input_file, "\n")
  
  # Read the input file without headers
  data <- tryCatch({
    read.table(input_file, header = FALSE, sep = "\t", skip = 1, check.names = FALSE)
  }, error = function(e) {
    cat("Error reading input file:", e$message, "\n")
    return(NULL)
  })

  if (is.null(data)) {
    stop("Failed to read input file.")
  }
  
  if (nrow(data) < 1) {
    cat("No data found in input file.")
    return()
  }
  
  # Check if the number of columns matches the expected column names
  if (length(col_names) != ncol(data)) {
    cat("Unexpected number of columns in input file:", ncol(data), "\n")
    stop("Unexpected number of columns in input file.")
  }
  
  # Assign the column names
  colnames(data) <- col_names
  
  # Debug: Print the first few rows and column names of the input file
  cat("First few rows of input data:\n")
  print(head(data))
  
  cat("Column names of input data:\n")
  print(colnames(data))
  
  # Ensure the data has the correct columns
  if (!all(required_columns %in% colnames(data))) {
    cat("Required columns not found in data. Columns present:", colnames(data), "\n")
    stop("Input file does not have the required columns.")
  }
  
  # Convert p1 to numeric and check for NAs if existing
  if ("p1" %in% colnames(data)) {
    data$p1 <- as.numeric(as.character(data$p1))
    if (any(is.na(data$p1))) {
      cat("NA values found in p1 column after conversion to numeric.\n")
      stop("NA values found in p1 column after conversion to numeric.")
    }
  }
  
  # Read the normalization file
  norm_data <- tryCatch({
    read.csv(norm_file, header = TRUE)
  }, error = function(e) {
    cat("Error reading normalization file:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(norm_data)) {
    stop("Failed to read normalization file.")
  }
  
  # Debug: Print column names and first few rows of the normalization data
  cat("Column names in normalization data:", colnames(norm_data), "\n")
  cat("First few rows of normalization data:\n")
  print(head(norm_data))
  
  # Convert bin files to numeric
  norm_data$bin <- as.numeric(as.character(norm_data$bin))
  norm_data$mean <- as.numeric(as.character(norm_data$mean))
  norm_data$std <- as.numeric(as.character(norm_data$std))
  
  # Function to find the closest bin
  find_closest_bin <- function(daf, norm_data) {
    idx <- which.min(abs(norm_data$bin - daf))
    return(norm_data[idx, ])
  }
  
  # Normalize the statistic
  data[[paste0("norm_", value_col)]] <- apply(data, 1, function(row) {
    daf_value <- as.numeric(row[["p1"]])  # Update this line to use "p1" instead of "daf"
    bin_info <- find_closest_bin(daf_value, norm_data)
    mean <- bin_info$mean
    std <- bin_info$std
    if (is.na(mean) || is.na(std) || std == 0) {
      return(NA)
    }
    norm_value <- (as.numeric(row[[value_col]]) - mean) / std
    return(round(norm_value, 4))
  })
  
  # Save the output
  write.table(data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# Helper function to convert value column name to lowercase for normalization file
normalize_value_column_name <- function(value_col) {
  return(tolower(value_col))
}

# Define the base path and file patterns
base_path <- "/home/tx56/palmer_scratch/deepsweep_empirical"
patterns <- list(
  list(pattern = "1kgp_xpehh", file_ext = ".joint.xpehh", value_col = "XPEHH", col_names = c("Index", "ID", "Freq", "iHH_A1", "iHH_B1", "iHH_P1", "XPEHH", "Std_XPEHH"), required_columns = c("Freq", "XPEHH")),
  list(pattern = "1kgp_ihs", file_ext = ".ihs", value_col = "iHS", col_names = c("Index", "ID", "Freq", "iHH_0", "iHH_1", "iHS", "Std_iHS"), required_columns = c("Freq", "iHS")),
  list(pattern = "1kgp_nsl", file_ext = ".nsl.out", value_col = "nsl", col_names = c("rs", "pos", "daf", "sl1", "sl0", "nsl"), required_columns = c("daf", "nsl")),
  list(pattern = "1kgp_ihh12", file_ext = ".ihh12.out", value_col = "ihh12", col_names = c("id", "pos", "p1", "ihh12"), required_columns = c("p1", "ihh12"))
)

# Process each pattern
for (pattern_info in patterns) {
  pattern_dir <- file.path(base_path, pattern_info$pattern)
  files <- list.files(pattern_dir, pattern = paste0(".*", pattern_info$file_ext), full.names = TRUE)
  
  for (file in files) {
    file_name <- basename(file)
    parts <- strsplit(file_name, "\\.")[[1]]
    pop_name <- parts[1]
    chr_num <- parts[2]
    
    # Construct the normalization file path with the value column in lowercase
    norm_file <- file.path(base_path, "bins", paste0(pop_name, "_", tolower(pattern_info$value_col), "_bin.csv"))
    output_file <- file.path(output_dir, paste0("norm_", pattern_info$value_col, "_", pop_name, "_", chr_num, ".tsv"))
    
    # Check if the output file already exists
    if (file.exists(output_file)) {
      cat("Output file already exists, skipping:", output_file, "\n")
      next
    }
    
    normalize_data(file, norm_file, output_file, pattern_info$value_col, pattern_info$col_names, pattern_info$required_columns)
  }
}

cat("Normalization completed.\n")
