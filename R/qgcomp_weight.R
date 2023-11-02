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
qgcomp_extract_weights <- function(folder_location) {
  # get list of files that start with "nb"
  nb_files <- list.files(path = folder_location, pattern = "^nb.*\\.rda$", full.names = TRUE)

  # create empty list to store results
  results_list <- list()

  # loop through each nb file
  for (file in nb_files) {
    # load data from file
    load(file)

    # extract neg and pos weights
    neg <- as.data.frame(nb$neg.weights) |>
      rename(weight = `nb$neg.weights`) |>
      rownames_to_column("chemical") |>
      mutate(direction="neg")

    pos <- as.data.frame(nb$pos.weights) |>
      rename(weight = `nb$pos.weights`) |>
      rownames_to_column("chemical") |>
      mutate(direction="pos")

    weights_df <- bind_rows(neg, pos)

    # extract file name without path and .rda extension using regular expression
    file_name <- sub(".*/(.*)\\.rda", "\\1", file)

    # add file name as a variable to weights_df
    weights_df$file_name <- file_name

    # append weights_df to results_list
    results_list[[file_name]] <- weights_df
  }

  # combine all results into a single dataframe
  results_df <- do.call(rbind, results_list)

  return(results_df)
}
