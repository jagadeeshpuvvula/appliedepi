#' mwas loop - perform linear regression between each feature and exposure metric provided in the input
#' metabolome data will be transposed in the function
#' input metabolome features and exposures are in different dataframes
#' runtime 17,000 * 70 * 17 matrix takes around 10 mins
#' @param feature_table provide a metabolome dataframe where each observation is a metabolome feature and each column is a subject
#' @param exp_cov_data each row represent a subject and contains exposure and covariate variables as columns
#' @param output_folder provide the location to the folder to save the output file from this function
#' @param mwas_file_name name of the CSV file with the linear regression results that will be save in output_folder
#' @param fdr_cutoff specify the FDR cutoff value here
#' @param exposures provide the variable names that represent the exposures
#' @param covar provide list of covariates that need to be adjusted in the model
#' @return This function will return two CSV files in the output_folder location provided by user 1st file contain linear regression results and 2nd with FDR threshold values
#' @export
#' @examples
#' \donttest{
#'mwas_lm_loop(feature_table = mom_ft,
#'             exp_cov_data = mom_dat,
#'             exposures = names(mom_dat)[4:20],
#'             covar = c("bmi", "mom_age", "pgtob", "mom_edu", "mom_race"),
#'             output_folder = "~/Documents/phth_phe_MWAS/result/lm_mwas",
#'             mwas_file_name = "mom_mwas.csv", #computes lm by each exposure but compiles all results in a single csv
#'             fdr_cutoff = 0.2)
#' }

mwas_lm_loop <- function(feature_table, exp_cov_data, output_folder, mwas_file_name, cutoff_file_name, exposures, covar, fdr_cutoff) {
  result_list <- list() # create an empty list to store results

  for (variable in exposures) {
    current_results <- data.frame(matrix(nrow = nrow(feature_table), ncol = 4)) # create a new data frame for each variable
    names(current_results) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    for (i in 1:nrow(feature_table)) {
      metabolite <- unlist(feature_table[i, -c(1)]) # store the ith metabolite in a new vector
      formula_str <- paste("metabolite ~", variable, "+", paste(covar, collapse = " + "))
      current_results[i, ] <- summary(lm(formula_str, data = cbind(metabolite, exp_cov_data)))$coefficients[2, c(1:4)] # store regression results for the ith metabolite
    }
    current_results$FDR <- p.adjust(current_results[[4]], method = "fdr") # fdr adjustment
    current_results <- cbind(feature_table[c(1)], current_results) # add mz and time
    current_results$Variable <- variable
    result_list[[variable]] <- current_results # add current results to list
  }

  PAH_MWAS_result <- do.call(rbind, result_list) # combine all results into one data frame

  # Calculate beta_dir variable
  PAH_MWAS_result$beta_dir <- case_when(
    PAH_MWAS_result$FDR < fdr_cutoff & PAH_MWAS_result$Estimate > 0 ~ "positive-significant",
    PAH_MWAS_result$FDR < fdr_cutoff & PAH_MWAS_result$Estimate < 0 ~ "negative-significant",
    PAH_MWAS_result$FDR > fdr_cutoff & PAH_MWAS_result$Estimate > 0 ~ "positive-non_significant",
    PAH_MWAS_result$FDR > fdr_cutoff & PAH_MWAS_result$Estimate < 0 ~ "negative-non_significant", TRUE ~ NA_character_)

  # Save final result to a file
  output_file <- file.path(output_folder, mwas_file_name)
  write.csv(PAH_MWAS_result, file = output_file, row.names = FALSE)

  return(PAH_MWAS_result)
}
