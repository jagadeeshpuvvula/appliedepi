#' Extract bkmr PIPs from all RDA files - prep for Shiny object
#'
#' @param folder_location input a folder location that contain bkmr output in RDA format
#'
#' @return This function will return a dataframe that contain information from file names and PIPs
#' @export
#'
#' @import bkmr
#'
#' @example
#' \donttest{
#' pip_df<- extract_pips("E:/BBK17/pj/data_2023apr/results/bkmr/lasso_15chem")
#' }
#'
extract_pips <- function(folder_location) {
  rda_files <- list.files(folder_location, pattern = "\\.rda$", full.names = TRUE)
  pips_list <- lapply(rda_files, function(file) {
    load(file)
    pip <- as.data.frame(ExtractPIPs(bkmr))
    pip$file_name <- sub(".*/(.*)\\.rda", "\\1", file)
    pip
  })
  pips_df <- do.call(rbind, pips_list)
  return(pips_df)
}
