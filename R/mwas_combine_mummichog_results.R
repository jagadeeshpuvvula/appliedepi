#' Combines mummuchog output files into a single dataframe
#'
#' @param folder_location Specify folder location that contain all the mummichog results
#'
#' @return returns a dataframe with all results into the global environment
#' @export
#'
#' @examples
#' \donttest{
#'enrich_data <- combine_csv_files("~/Documents/phth_phe_MWAS/result/enrich_imp")
#' }
#'
#'
combine_csv_files <- function(folder_location) {
  # Get list of CSV files in the folder
  file_list <- list.files(path = folder_location, pattern = "\\.csv$", full.names = TRUE)

  # Initialize an empty list to store data frames
  data_frames <- list()

  # Read and row bind each CSV file
  for (file in file_list) {
    # Extract chemical and mom_baby from file name
    file_name <- basename(file)
    file_parts <- strsplit(file_name, "_")[[1]]
    chemical <- file_parts[1]
    mom_baby <- gsub(".csv", "", file_parts[2], fixed = TRUE)

    # Read CSV file
    csv_data <- read.csv(file)

    # Add chemical and mom_baby as variables
    csv_data$chemical <- chemical
    csv_data$mom_baby <- mom_baby

    # Append the data frame to the list
    data_frames[[length(data_frames) + 1]] <- csv_data
  }

  # Combine all data frames into a single data frame
  combined_df <- do.call(rbind, data_frames)

  return(combined_df)
}
