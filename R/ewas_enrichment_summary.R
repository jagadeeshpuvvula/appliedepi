#' ewas enrichemnt summary
#'
#' @param folder_path folder location that contain the enriched results in csv format
#' @param output_folder folder to save results
#'
#' @return
#' @export
#'
#' @example
#' \donttest{
enrichment_summary(folder_path = "~/Documents/methylation_pilot/res_enriched_cell_cnt",
                   output_folder = "~/Documents/methylation_pilot/res_enriched_cell_cnt")
#' }
enrichment_summary <- function(folder_path, output_folder) {
  # Get a list of CSV files in the folder
  csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

  # Function to read a CSV file, select top 10 smallest P.DE values, and add the file name
  read_and_select_top_n <- function(file_path) {
    file_name <- tools::file_path_sans_ext(basename(file_path))

    # Split the filename based on underscores and dots
    name_parts <- unlist(strsplit(file_name, "[._]"))

    # Extract chemical and tissue
    chemical <- paste(name_parts[-length(name_parts)], collapse = "_")
    tissue <- name_parts[length(name_parts)]

    # Read CSV and add "pathway_id" as row names
    df <- read.csv(file_path, header = TRUE)
    df <- df[order(df$P.DE), ][1:10, ]
    df$chemical <- chemical
    df$tissue <- tissue

    return(df)
  }

  # Read and process each CSV file, and bind the results into a single data frame
  results <- lapply(csv_files, read_and_select_top_n)

  # Combine the results into a single data frame
  results_df <- do.call(rbind, results)

  # Export the results to the output folder
  output_file <- file.path(output_folder, "enrichment_summary.csv")
  write.csv(results_df, file = output_file, row.names = TRUE)

  return(output_file)
}
