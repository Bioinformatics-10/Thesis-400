# 1. Get the path (e.g., from your file explorer)
my_folder_path <- "A:/thesis/network_building"
# OR 
# my_folder_path <- "/Users/MyName/Documents/Gene_Sets" 

# 2. Set the working directory
setwd(my_folder_path)

# 3. Verify it worked
getwd()


gmt_to_txt <- function(gmt_file, output_file) {
  # Check if the GMT file exists before attempting to read
  if (!file.exists(gmt_file)) {
    stop(paste("Error: GMT file not found at path:", gmt_file))
  }
  
  # Read the GMT file line by line (Base R)
  gmt_lines <- readLines(gmt_file)
  
  # Process each line to extract genes
  all_genes <- unique(unlist(lapply(gmt_lines, function(x) {
    # Split the line by tab character (Base R)
    parts <- strsplit(x, "\t")[[1]]
    # Genes start from the 3rd element (index 3) after the set name and description
    genes <- parts[-c(1, 2)]
    return(genes)
  })))
  
  # Write the unique genes to a plain text file (Base R)
  write.table(all_genes, file = output_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Display a completion message (Base R)
  message("Finished: All unique genes extracted and written to ", output_file)
}

# Example Usage:
# NOTE: Ensure the file 'c7.all.v2024.1.Hs.symbols.gmt' is in your R working directory.
gmt_to_txt("A:/thesis/network_building/c7.all.v2025.1.Hs.symbols.gmt", "C7_genes.txt")

#immport_genes <- read.table("A:/thesis/network_building/all_gene_lists.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]

# This script merges two gene lists from two separate files, ensuring 
# all gene names are unique in the final output.

# --- Step 1: Define File Paths ---
# Use the full, absolute path for the first file. 
# IMPORTANT: In R, use forward slashes (/) for all paths, even on Windows.
path_to_all_gene_lists <- "A:/thesis/network_building/all_gene_lists.txt"

# For the second file, 'C7_genes.txt', we assume it is in the current R working directory (getwd()).
# If it is not, replace the simple filename with its full absolute path.
path_to_c7_genes <- "C7_genes.txt" 

output_file_name <- "Merged_Immune_Genes.txt"

# --- Step 2: Load the First Gene List (all_gene_lists.txt) ---
tryCatch({
  # read.table is used assuming one gene symbol per line.
  # We use sep = "\n" to ensure it reads the entire line as a single field
  # in case there are spaces or odd characters, and then extract the first column ([[1]])
  all_gene_lists <- read.table(path_to_all_gene_lists, 
                               header = FALSE, 
                               stringsAsFactors = FALSE, 
                               sep = "\n")[[1]]
  message("Successfully loaded: ", path_to_all_gene_lists)
}, error = function(e) {
  stop(paste("Error loading", path_to_all_gene_lists, ":", conditionMessage(e)))
})


# --- Step 3: Load the Second Gene List (C7_genes.txt) ---
tryCatch({
  c7_genes <- read.table(path_to_c7_genes, 
                         header = FALSE, 
                         stringsAsFactors = FALSE, 
                         sep = "\n")[[1]]
  message("Successfully loaded: ", path_to_c7_genes)
}, error = function(e) {
  # If this file is missing, it suggests the gmt_processor.R script may not have run, 
  # or R's working directory is incorrect.
  stop(paste("Error loading", path_to_c7_genes, ":", conditionMessage(e),
             "\nNote: This file is expected to be in your current working directory:", getwd()))
})


# --- Step 4: Merge and Find Unique Genes ---
# Combine the two character vectors and find the unique set of genes (Base R)
merged_immune_genes <- unique(c(all_gene_lists, c7_genes))

# --- Step 5: Save the Final List ---
# The new file will be saved in the current R working directory (getwd())
write.table(merged_immune_genes, output_file_name,
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

# --- Step 6: Output Summary ---
message("Finished merging gene lists.")
message("Total unique genes found: ", length(merged_immune_genes))
message("Output file saved to current directory: ", file.path(getwd(), output_file_name))

# Optional: Display the total count
length(merged_immune_genes)