#' This function will bind rows across the csv files in a given target folder
#' @param folder_path User need to prove the folder location that contain csv files that end with mwas
#' these csv files with mwas as a prefix will be resulted from the mwas_lm_loop function above
#' Especially helpful if the results are in multiple files
#' @export
#' @return this function will return a dataframe within the global environment
#' @examples
#' \donttest{
#' data_combined <- bind_mwas_res_from_folder(paste0(result,"lm_mwas"))
#' }
#'
#'

bind_mwas_res_from_folder <- function(folder_path) {
  # Get a list of all CSV files in the folder ending with "mwas"
  csv_files <- list.files(path = folder_path, pattern = "mwas\\.csv$", full.names = TRUE)

  # Create an empty list to store data frames from individual CSV files
  df_list <- vector("list", length = length(csv_files))

  # Read each CSV file and store it in the list
  for (i in seq_along(csv_files)) {
    df_list[[i]] <- read.csv(csv_files[i])
  }

  # Combine the data frames from the list row-wise
  data_combined <- do.call(rbind, df_list)

  return(data_combined)
}
