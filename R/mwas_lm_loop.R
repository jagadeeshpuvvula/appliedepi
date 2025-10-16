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
#'mwas_lm_loop(feature_table = mom_ft, #transposed feature table - each observation is a metabolite
#'             exp_cov_data = mom_dat,
#'             exposures = names(mom_dat)[4:20],
#'             covar = c("bmi", "mom_age", "pgtob", "mom_edu", "mom_race"),
#'             output_folder = "~/Documents/phth_phe_MWAS/result/lm_mwas",
#'             mwas_file_name = "mom_mwas.csv", #compiles all results in a single csv
#'             fdr_cutoff = 0.2,
#'             sex_var = "gender", analyze_by_sex = TRUE, participant_id_col = "participant_id")
#' }

mwas_lm_loop <- function(feature_table, exp_cov_data, output_folder, mwas_file_name, cutoff_file_name, exposures, covar, fdr_cutoff, sex_var = NULL, analyze_by_sex = FALSE, participant_id_col = "participant_id") {
  
  # ===== QC STEP 1: Validate input parameters =====
  if (analyze_by_sex && is.null(sex_var)) {
    stop("sex_var must be provided when analyze_by_sex is TRUE")
  }
  
  if (!participant_id_col %in% names(exp_cov_data)) {
    stop(paste("participant_id column '", participant_id_col, "' not found in exp_cov_data", sep = ""))
  }
  
  # ===== QC STEP 2: Validate dimensions =====
  message(paste("QC: feature_table dimensions:", nrow(feature_table), "metabolites x", ncol(feature_table) - 1, "samples"))
  message(paste("QC: exp_cov_data dimensions:", nrow(exp_cov_data), "participants x", ncol(exp_cov_data), "variables"))
  
  # ===== QC STEP 3: Extract and match participant IDs =====
  feature_table_ids <- names(feature_table)[2:ncol(feature_table)]
  exp_cov_ids <- exp_cov_data[[participant_id_col]]
  
  message(paste("QC: feature_table has", length(feature_table_ids), "participants"))
  message(paste("QC: exp_cov_data has", length(exp_cov_ids), "participants"))
  
  # Check for overlapping IDs
  matching_ids <- intersect(feature_table_ids, exp_cov_ids)
  message(paste("QC: Found", length(matching_ids), "matching participants"))
  
  if (length(matching_ids) == 0) {
    stop("No matching participant IDs found between feature_table and exp_cov_data")
  }
  
  if (length(matching_ids) < length(exp_cov_ids) * 0.8) {
    warning(paste("Warning: Less than 80% of participants matched. Matching rate:", 
                  round(length(matching_ids) / length(exp_cov_ids) * 100, 1), "%"))
  }
  
  # ===== QC STEP 4: Reorder and subset data to match =====
  exp_cov_data <- exp_cov_data[match(matching_ids, exp_cov_ids), ]
  feature_table <- feature_table[, c(1, match(matching_ids, feature_table_ids) + 1)]
  
  if (!all(names(feature_table)[2:ncol(feature_table)] == exp_cov_data[[participant_id_col]])) {
    stop("Failed to properly match participant IDs between datasets")
  }
  
  message(paste("QC: Successfully matched and reordered", nrow(exp_cov_data), "participants"))
  
  # ===== QC STEP 5: Check for missing values =====
  missing_exp_cov <- colSums(is.na(exp_cov_data))
  if (any(missing_exp_cov > 0)) {
    message("QC: Missing values in exp_cov_data:")
    print(missing_exp_cov[missing_exp_cov > 0])
  }
  
  missing_features <- rowSums(is.na(feature_table[, -1]))
  if (any(missing_features > 0)) {
    message(paste("QC: Found", sum(missing_features > 0), "metabolites with missing values"))
  }
  
  # ===== QC STEP 6: Validate exposures and covariates exist =====
  missing_exposures <- exposures[!exposures %in% names(exp_cov_data)]
  if (length(missing_exposures) > 0) {
    stop(paste("Missing exposures:", paste(missing_exposures, collapse = ", ")))
  }
  
  missing_covar <- covar[!covar %in% names(exp_cov_data)]
  if (length(missing_covar) > 0) {
    stop(paste("Missing covariates:", paste(missing_covar, collapse = ", ")))
  }
  
  message(paste("QC: All", length(exposures), "exposures and", length(covar), "covariates found"))
  
  # ===== QC STEP 7: Validate sex_var if needed =====
  if (analyze_by_sex) {
    if (!sex_var %in% names(exp_cov_data)) {
      stop(paste("sex_var '", sex_var, "' not found in exp_cov_data", sep = ""))
    }
    sex_summary <- table(exp_cov_data[[sex_var]], useNA = "ifany")
    message("QC: Sex variable distribution:")
    print(sex_summary)
  }
  
  result_list <- list()
  
  # Helper function to run analysis
  run_analysis <- function(current_feature_table, current_exp_cov_data, current_covar, analysis_name, sex_group_label) {
    message(paste("=== Analyzing:", analysis_name, "(n =", nrow(current_exp_cov_data), ") ==="))
    
    for (variable in exposures) {
      current_results <- data.frame(matrix(nrow = nrow(current_feature_table), ncol = 4))
      names(current_results) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      
      for (i in 1:nrow(current_feature_table)) {
        metabolite <- unlist(current_feature_table[i, -c(1)])
        
        # Create formula
        covar_str <- if (length(current_covar) > 0) paste("+", paste(current_covar, collapse = " + ")) else ""
        formula_str <- as.formula(paste("metabolite ~", variable, covar_str))
        
        # Create model data
        model_data <- cbind(metabolite = metabolite, current_exp_cov_data)
        
        tryCatch({
          model <- lm(formula_str, data = model_data)
          current_results[i, ] <- summary(model)$coefficients[2, c(1:4)]
        }, error = function(e) {
          current_results[i, ] <<- NA
          warning(paste("Model fitting failed for metabolite", i, "and variable", variable))
        })
      }
      
      current_results$FDR <- p.adjust(current_results[[4]], method = "fdr")
      current_results <- cbind(current_feature_table[c(1)], current_results)
      current_results$Variable <- variable
      
      if (analyze_by_sex) {
        current_results$Sex <- sex_group_label
      }
      
      result_key <- paste(variable, sex_group_label, sep = "_")
      result_list[[result_key]] <<- current_results
    }
  }
  
  # ===== STEP 1: OVERALL ANALYSIS (all observations with sex_var as covariate only if analyzing by sex) =====
  message("\n>>> STEP 1: OVERALL ANALYSIS <<<")
  
  # When performing sex-stratified analysis, include sex_var as a covariate in the overall model to adjust for sex effects
  if (analyze_by_sex) {
    overall_covar <- unique(c(covar, sex_var))
    overall_sex_label <- "overall"
  } else {
    overall_covar <- covar
    overall_sex_label <- "all"
  }
  
  run_analysis(
    current_feature_table = feature_table,
    current_exp_cov_data = exp_cov_data,
    current_covar = overall_covar,
    analysis_name = "Overall (all observations)",
    sex_group_label = overall_sex_label
  )
  
  # ===== STEP 2: SEX-STRATIFIED ANALYSIS (if analyze_by_sex = TRUE) =====
  if (analyze_by_sex) {
    message("\n>>> STEP 2: SEX-STRATIFIED ANALYSIS <<<")
    sex_levels <- unique(exp_cov_data[[sex_var]])
    sex_levels <- sex_levels[!is.na(sex_levels)]
    sex_levels <- sort(sex_levels)
    
    for (sex_group in sex_levels) {
      group_indices <- which(exp_cov_data[[sex_var]] == sex_group)
      current_exp_cov_data <- exp_cov_data[group_indices, ]
      current_feature_table <- feature_table[, c(1, group_indices + 1)]
      current_covar <- covar[covar != sex_var]  # Remove sex_var for stratified analysis
      
      run_analysis(
        current_feature_table = current_feature_table,
        current_exp_cov_data = current_exp_cov_data,
        current_covar = current_covar,
        analysis_name = paste(sex_group, "only"),
        sex_group_label = sex_group
      )
    }
  }
  
  MWAS_result <- do.call(rbind, result_list)
  rownames(MWAS_result) <- NULL
  
  # Calculate beta_dir variable
  MWAS_result$beta_dir <- dplyr::case_when(
    MWAS_result$FDR < fdr_cutoff & MWAS_result$Estimate > 0 ~ "positive-significant",
    MWAS_result$FDR < fdr_cutoff & MWAS_result$Estimate < 0 ~ "negative-significant",
    MWAS_result$FDR >= fdr_cutoff & MWAS_result$Estimate > 0 ~ "positive-non_significant",
    MWAS_result$FDR >= fdr_cutoff & MWAS_result$Estimate < 0 ~ "negative-non_significant",
    TRUE ~ NA_character_
  )
  
  # Save final result to a file
  output_file <- file.path(output_folder, mwas_file_name)
  write.csv(MWAS_result, file = output_file, row.names = FALSE)
  
  message("\n>>> ANALYSIS COMPLETE <<<")
  message(paste("Results saved to:", output_file))
  
  return(MWAS_result)
}
