#' Summarize results from EWAS analysis
#'
#' @param folder_path provide location of a folder that contains EWAS results in csv format
#'
#' @return returns a summary of the associations by exposure-methylation data combinations
#' @export
#'
#' @example
#' \donttest{
#' ewas_summary("~/Documents/methylation_pilot/res_cell_cnt")
#' }
ewas_summary <- function(folder_path) {
  # Get a list of CSV files in the specified folder
  files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

  # Initialize a parallel backend with the number of available CPU cores
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  # Perform file processing in parallel
  result_list <- foreach(file_path = files) %dopar% {
    # Read the file
    data <- read.csv(file_path)

    # Filter observations with FDR < 0.05
    filtered_data <- data[data$FDR < 0.05, ]

    # Get the count of observations with FDR < 0.05
    count_fdr_lt_0.05 <- nrow(filtered_data)

    # Get the list of CPG.Labels with FDR < 0.05 as a comma-separated string
    cpg_labels_lt_0.05 <- paste(filtered_data$CPG.Labels, collapse = ", ")

    # Extract file name from the path
    file_name <- basename(file_path)

    # Return the result for the current file
    return(list(
      File_Name = file_name,
      Num_Observations_FDR_lt_0.05 = count_fdr_lt_0.05,
      CPG_Labels_FDR_lt_0.05 = cpg_labels_lt_0.05
    ))
  }

  # Stop the parallel backend
  stopCluster(cl)

  # Convert the list of results to a data frame
  result_df <- do.call(rbind.data.frame, result_list)

  # Return the resulting dataframe
  return(result_df)
}
