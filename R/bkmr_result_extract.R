#' Extract results from the bkmr RDA files
#'
#' @param folder_location input folder location that contain bkmr output in RDA format
#' @param select_iterations select iterations among the total number of iterations
#' @param estimate_sequence input sequence and intervals to display overall estimates.
#'
#' @return this function will return a list of two objects the first is PIPs and second is overall risk estimates
#' @export
#'
#' @import bkmr
#'
#' @example
#' #' \donttest{
#'results <- extract_data_from_folder("E:/BBK17/pj/data_2023apr/results/bkmr/lasso_15chem",
#'                                    select_iterations= seq(25000, 50000, by = 50),
#'                                    estimate_sequence= seq(0.10, 0.90, by = 0.05))
#' }
extract_data_from_folder <- function(folder_location, select_iterations, estimate_sequence) {

  # Get a list of all .rda files in the folder
  rda_files <- list.files(folder_location, pattern = ".rda$")

  # Initialize empty data frames for pips and overall risk
  pips_df <- data.frame()
  overall_risk_df <- data.frame()

  # Loop through each .rda file in the folder
  for (rda_file in rda_files) {
    # Load the bkmr object from the .rda file
    load(file.path(folder_location, rda_file))

    # Split the file name into its component parts
    file_parts <- strsplit(rda_file, "_")[[1]]
    model_name <- file_parts[1]
    data_name <- file_parts[2]
    outcome_name <- file_parts[3]
    covariate_names <- paste(file_parts[4:(length(file_parts)-1)], collapse = "_")
    measure_name <- gsub(".rda", "", file_parts[length(file_parts)])

    # Extract the PIPs and add them to the pips_df data frame
    pips <- ExtractPIPs(bkmr) %>% as.data.frame()
    pips$model_name <- model_name
    pips$data_name <- data_name
    pips$outcome_name <- outcome_name
    pips$covariate_names <- covariate_names
    pips$measure_name <- measure_name
    pips_df <- bind_rows(pips_df, pips)

    # Extract the overall risk and add it to the overall_risk_df data frame
    sel <- select_iterations
    overall_risk <- OverallRiskSummaries(fit = bkmr, qs = estimate_sequence,
                                         q.fixed = 0.5, method = "approx", sel = sel)
    overall_risk$model_name <- model_name
    overall_risk$data_name <- data_name
    overall_risk$outcome_name <- outcome_name
    overall_risk$covariate_names <- covariate_names
    overall_risk$measure_name <- measure_name
    overall_risk_df <- bind_rows(overall_risk_df, overall_risk)
  }

  # Return the pips and overall risk data frames
  return(list(pips_df = pips_df, overall_risk_df = overall_risk_df))
}
