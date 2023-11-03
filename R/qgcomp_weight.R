#' Extract chemical weights contributed to the joint association using non bootstrap quantile g-computation.
#'
#' @param folder_location input folder location that contain the RDA output files from the non bootstrap quantile g computation. this function will extract from files that contain nb in the filename.
#'
#' @return returns a dataframe that will contain the RDA fileaname, chemical name and corresponding weight
#' @export
#'
#' @import tidyverse
#'
#' @example
#' \donttest{
#'qgcomp_weights <- qgcomp_extract_weights("E:/BBK17/pj/data_2023apr/results/qgcomp/qgcomp_pcs_wo_wisc")
#' }
extract_weights <- function(folder_location) {
  # get list of files that start with "nb"
  nb_files <- list.files(path = folder_location, pattern = "^nb.*\\.rda$", full.names = TRUE)
  
  # create an empty list to store results
  results_list <- list()
  
  # loop through each nb file
  for (file in nb_files) {
    # load data from file
    load(file)
    
    # extract neg and pos weights
    neg <- as.data.frame(nb$neg.weights) |> 
      rename(weight = `nb$neg.weights`) |> 
      rownames_to_column("chemical") |> 
      mutate(direction = "neg")
    
    pos <- as.data.frame(nb$pos.weights) |> 
      rename(weight = `nb$pos.weights`) |> 
      rownames_to_column("chemical") |> 
      mutate(direction = "pos")
    
    weights_df <- bind_rows(neg, pos)
    
    # extract file name without path and .rda extension using regular expression
    file_name <- sub(".*/(.*)\\.rda", "\\1", file)
    
    # Add "file_name" as a variable to weights_df
    weights_df$file_name <- file_name
    
    # Extract psi values and add them as rows
    if (is.list(nb$psi)) {
      psi1 <- nb$psi$psi1
    } else {
      psi1 <- NA
    }
    
    psi_df <- data.frame(chemical = "neg.psi", weight = nb$neg.psi, direction = "neg", file_name = file_name) |>
      rbind(data.frame(chemical = "pos.psi", weight = nb$pos.psi, direction = "pos", file_name = file_name)) |>
      rbind(data.frame(chemical = "psi1", weight = nb$psi, direction = "psi1", file_name = file_name))
    
    # Append psi_df to weights_df
    combined_df <- bind_rows(weights_df, psi_df)
    
    # Append combined_df to results_list
    results_list[[file_name]] <- combined_df
  }
  
  # Combine all results into a single dataframe
  results_df <- do.call(rbind, results_list)
  
  return(results_df)
}
