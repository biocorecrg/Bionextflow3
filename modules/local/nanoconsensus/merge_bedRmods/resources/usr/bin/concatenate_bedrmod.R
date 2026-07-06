#!/usr/bin/env Rscript

# Concatenate bedrmod files, merging modification_names and appending data rows
# Usage: Rscript concatenate_bedrmod.R file1.bedrmod file2.bedrmod ... output.bedrmod

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript concatenate_bedrmod.R file1.bedrmod [file2.bedrmod ...] output.bedrmod\n")
  cat("At least 1 input file and 1 output file must be provided.\n")
  quit(status = 1)
}

# Separate input files from output file
output_file <- args[length(args)]
input_files <- args[-length(args)]

# Find first non-empty input file to use as template
first_file <- NULL
first_lines <- NULL

for (f in input_files) {
  lines <- readLines(f)
  if (length(lines) > 0 && nzchar(lines[1])) {
    first_file <- f
    first_lines <- lines
    break
  }
}

if (is.null(first_file)) {
  message("All input files are empty — nothing to merge, exiting.")
  quit(save = "no", status = 0)
}

# Determine input file (bedRmod/bed):
is_track_format <- grepl("^track name=", first_lines[1])

if (is_track_format) {
  # Simple track format: header is just the first line, no modification merging needed
  header_lines <- ""
  output_extension <- "bed"
  column_header <- first_lines[1]

} else {

    # Find where header ends (first line not starting with #)
    header_end <- which(!grepl("^#", first_lines))[1] - 1
    header_lines <- first_lines[1:header_end]
    column_header <- first_lines[header_end + 1]

    # Extract all modification names from all files
    all_modifications <- c()

    for (file in input_files) {
      lines <- readLines(file)
      mod_line <- lines[grepl("^#modification_names=", lines)]

      if (length(mod_line) > 0) {
        # Extract the modification names value
        mod_names <- sub("^#modification_names=", "", mod_line)
        # Split by comma and add to all_modifications
        mods <- strsplit(mod_names, ",")[[1]]
        all_modifications <- c(all_modifications, mods)
      }
    }

    # Get unique modification names and sort them
    unique_modifications <- unique(all_modifications)
    unique_modifications <- sort(unique_modifications)
    combined_mod_names <- paste(unique_modifications, collapse = ",")

    # Update the modification_names line in header
    header_lines <- gsub("^#modification_names=.*",
                         paste0("#modification_names=", combined_mod_names),
                         header_lines)

    output_extension <- "bedrmod"

}

# Collect all data rows (non-header lines) from all files
all_data_rows <- c()

for (file in input_files) {
  lines <- readLines(file)
  # Skip all header lines (starting with #) and the column header line
  data_rows <- lines[!grepl("^#", lines) & lines != column_header & lines != ""]
  all_data_rows <- c(all_data_rows, data_rows)
}

# Write output file
output_lines <- c(header_lines, column_header, all_data_rows)
writeLines(output_lines, sprintf("%s.%s", output_file, output_extension))

cat(sprintf("Successfully concatenated %d files into %s\n",
            length(input_files), output_file))
cat(sprintf("Total data rows: %d\n", length(all_data_rows)))
