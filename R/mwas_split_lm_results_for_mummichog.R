#' Split MWAS LM result files from a single dataframe to .txt files for PEA mummichog
#'
#' @param data dataframe name that contains all the MWAS LM results combined
#' @param variables specify exposure variable as variable and mom/baby variable as set
#' @param export_location specify the folder location to save all the .txt files ready for mummichog analysis
#'
#' @return returns .txt files at the specified export_location folder
#' @export
#'
#' @examples
#' \donttest{
#'split_and_export(data= dat,
#'                 variables = c("variable", "set"),
#'                 export_location = "~/Documents/phth_phe_MWAS/result/enrich_exp/")
#'
#' }
#'
split_and_export <- function(data, variables, export_location) {
  unique_combinations <- unique(data[, variables, drop = FALSE])

  for (i in 1:nrow(unique_combinations)) {
    filter_conditions <- paste(paste(names(unique_combinations), '=="', unique_combinations[i, ], '"', sep = ""), collapse = ' & ')
    subset_data <- data %>% filter(!!rlang::parse_expr(filter_conditions))

    # Construct filename based on variable and set names
    filename <- paste(unique_combinations[i, ], collapse = "_")
    filepath <- file.path(export_location, paste0(filename, ".txt"))

    write.table(subset_data, filepath, sep = "\t", row.names = FALSE)
  }
}
