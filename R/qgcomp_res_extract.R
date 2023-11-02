#' Extract results from the qgcomp output RDA files. Iterates though the RDA files and compiles results in tibble format
#'
#' @param results_folder input a folder location that contain qgcomp RDA output files
#'
#' @return results a dataframe that will contain bootstap or not, gender, cohort, outcome, estimate, confidence intervals
#' @export
#'
#' @example
#' \donttest{
#'results<- get_gcomp_estimates(results_folder = "E:/BBK17/pj/data_2023apr/results/qgcomp/qgcomp_pcs_wo_wisc")
#' }
#'
#'
#'
#'
get_gcomp_estimates<- function(results_folder) {
  # get the file names in the folder
  file_names <- list.files(results_folder)

  # create empty data frame to store the results
  results <- data.frame(file = character(),
                        estimate = numeric(),
                        lower_ci = numeric(),
                        upper_ci = numeric(),
                        stringsAsFactors = FALSE)

  #loops through all .rda files and extract estimates, CI & p-values
  for (file_name in file_names) {
    # load the object from the file
    load(file.path(results_folder, file_name))
    # extract the object name using grep or grepl
    object_name <- ls()[grep("^boot|^nb", ls())]
    # extract the required information using the summary function and the object name
    estimate <- get(object_name)$coef["psi1"]
    CI <- get(object_name)$ci
    lower_ci <- CI[1]
    upper_ci <- CI[2]
    p_value <- get(object_name)$pval[2]
    # split file name by "_" and remove ".rda"
    file_parts <- strsplit(gsub("\\.rda", "", file_name), "_")[[1]]
    # assign parts to variables
    part1 <- file_parts[1]
    part2 <- file_parts[2]
    part3 <- file_parts[3]
    part4 <- file_parts[5]
    # print out the values
    cat("File:", file_name, "\n")
    cat("  Estimate:", estimate, "\n")
    cat("  Lower CI:", lower_ci, "\n")
    cat("  Upper CI:", upper_ci, "\n")
    cat("  p-value:", p_value, "\n")
    cat("  Part 1:", part1, "\n")
    cat("  Part 2:", part2, "\n")
    cat("  Part 3:", part3, "\n")
    cat("  Part 4:", part4, "\n")
    # add the results to the data frame
    results <- rbind(results, data.frame(boot_strp = part1,
                                         gender = part2,
                                         cohort = part3,
                                         outcome = part4,
                                         estimate = estimate,
                                         lower_ci = lower_ci,
                                         upper_ci = upper_ci,
                                         p_value = p_value,
                                         stringsAsFactors = FALSE))
    # remove the object from the R environment to avoid conflicts with other objects (all objects saved as nb or boot)
    rm(list = object_name)
  }

  return(results)
}
