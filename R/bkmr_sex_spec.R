#' BKMR wrapper performing sex specific loops. Load Matrix and bindrcpp packages
#'
#' @param data input of a dataframe that contain all the variables mentioned in the wrapper
#' @param folder_path input folder location to save the output in RDA format
#' @param include_sex input TRUE or FASE to perform sex stratified analysis
#' @param ln_mixture_vars input list of chemicals that need to be considered as mixture
#' @param outcome_vars input list of outcome variables
#' @param covariate_vars input list of covariates
#' @param groups input grouping structure of the chemicals considered as mixture, these values should be integers
#' @param iter input a numeric value for number of iterations
#'
#' @return this function will return output from bkmr model in RDA format
#' @export
#'
#' @import bkmr
#' @import fields
#'
#' @example
#' \donttest{
#' chemicals_29<- c("Pb","Hg", "DMA", "DDE", "PBDE_47","PCB_118","PCB_138", "PCB_153","PCB_180",
#'                  "PFHxS", "PFOA", "PFOS", "BCEtP", "BDCPP", "DBuP", "DPhP", "TCS", "BPA",
#'                  "di_Ethyl_OP", "di_Methyl_OP", "MBP","MBZP","MCPP", "sigma_DEHP","MEP",
#'                  "MIBP","B_PB", "M_PB", "P_PB")
#'chem_group_29<- c(rep(1,times=3), rep(2,times=6), rep(3, times=3),
#'                  rep(4,times=4), rep(5,times=2), rep(6,times=2),
#'                  rep(7,times=6), rep(8,times=3))
#'outcomes<-c("wppsi_fsiq","wppsi_viq", "wppsi_piq")
#'covariates<- c("cotinine", "home_score_total", "mom_age", "cohort",
#'                "city", "sex", "race_bin", "mom_edu_cat", "parity_n" )
#'
#' bkmr_cohort(data=dat_mirec,
#'            folder_path = "E:/BBK17/pj/data_2023apr/results/bkmr/29_chem/mirec_sex",
#'            ln_mixture_vars= chemicals_29,
#'            groups=chem_group_29,
#'            outcome_vars= outcomes,
#'            covariate_vars= covariates,
#'            include_cohort = TRUE)
#' }
bkmr_sex <- function(data, folder_path, include_sex = TRUE,
                     ln_mixture_vars = NULL, outcome_vars = NULL, covariate_vars = NULL, groups,
                     iter = 10000, min_n = 50) {
  set.seed(2023)
  
  # Check and standardize sex variable if needed
  if (include_sex) {
    if (!"sex" %in% names(data)) {
      stop("'sex' variable not found in data")
    }
    
    # Print sex coding to help debug
    cat("Original sex variable distribution:\n")
    print(table(data$sex, useNA = "always"))
    
    # Standardize sex coding to match male/female text values
    data$sex_clean <- tolower(trimws(as.character(data$sex)))
    
    cat("\nCleaned sex distribution:\n")
    print(table(data$sex_clean, useNA = "always"))
  }
  
  # Check variable types before starting
  cat("\nChecking variable types...\n")
  all_vars <- c(ln_mixture_vars, outcome_vars, covariate_vars)
  for (var in all_vars) {
    var_class <- class(data[[var]])
    if (is.list(data[[var]]) && !is.data.frame(data[[var]])) {
      stop(paste0("Variable '", var, "' is a list. Please convert to numeric vector."))
    }
    if (!is.numeric(data[[var]]) && !is.integer(data[[var]])) {
      cat("Warning: Variable '", var, "' is ", var_class, ". Converting to numeric.\n")
      data[[var]] <- as.numeric(data[[var]])
    }
  }
  cat("Variable type check complete.\n\n")
  
  # Define sex levels based on include_sex parameter
  if (include_sex) {
    sex_levels <- c("male", "female")
  } else {
    sex_levels <- c("all")
  }
  
  # Loop through sex levels
  for (sex_level in sex_levels) {
    # Subset data based on sex level
    if (include_sex) {
      sex_data <- subset(data, sex_clean == sex_level)
      current_covariate_vars <- covariate_vars[!(covariate_vars %in% c("sex"))]
    } else {
      sex_data <- data
      current_covariate_vars <- covariate_vars
    }
    
    cat("\n", rep("=", 60), "\n", sep = "")
    cat("Analyzing:", sex_level, "\n")
    cat("Total observations before filtering:", nrow(sex_data), "\n")
    cat(rep("=", 60), "\n\n", sep = "")
    
    # Run bkmr for each outcome variable
    for (outcome_var in outcome_vars) {
      # Create subset with complete cases for this specific outcome
      vars_needed <- c(ln_mixture_vars, outcome_var, current_covariate_vars)
      analysis_data <- sex_data[, vars_needed, drop = FALSE]
      
      # Ensure all columns are numeric
      for (var in names(analysis_data)) {
        if (is.list(analysis_data[[var]]) && !is.data.frame(analysis_data[[var]])) {
          stop(paste0("Variable '", var, "' is a list in analysis_data. Cannot proceed."))
        }
        analysis_data[[var]] <- as.numeric(analysis_data[[var]])
      }
      
      analysis_data <- analysis_data[complete.cases(analysis_data), , drop = FALSE]
      
      # Report sample size
      cat("Outcome:", outcome_var, "| N =", nrow(analysis_data), "\n")
      
      # Check if we have enough data
      if (nrow(analysis_data) < min_n) {
        warning(paste0("Insufficient data for ", sex_level, " and outcome ", outcome_var, 
                      ". N = ", nrow(analysis_data), " (minimum required: ", min_n, "). Skipping this combination."))
        next
      }
      
      # Extract variables needed for bkmr
      ln_mixture <- as.matrix(analysis_data[, ln_mixture_vars, drop = FALSE])
      outcome <- as.numeric(analysis_data[[outcome_var]])
      
      if (length(current_covariate_vars) > 0) {
        covariates <- as.matrix(analysis_data[, current_covariate_vars, drop = FALSE])
      } else {
        covariates <- NULL
      }
      
      # Verify no list columns remain
      if (is.list(outcome) && !is.numeric(outcome)) {
        stop(paste0("Outcome '", outcome_var, "' could not be converted to numeric vector"))
      }
      
      # Create knots for bkmr
      n_knots <- max(10, round(0.15 * nrow(ln_mixture)))
      knots100 <- fields::cover.design(ln_mixture, nd = n_knots)$design
      
      cat("Running BKMR with", n_knots, "knots...\n")
      
      # Run BKMR
      bkmr <- kmbayes(y = outcome, Z = ln_mixture, X = covariates, iter = iter,
                      verbose = TRUE, varsel = TRUE,
                      groups = groups,
                      knots = knots100)
      
      # Save results
      save(bkmr, file = paste0(folder_path, "/bkmr_",
                               sex_level, "_",
                               outcome_var, ".rda"))
      
      cat("Saved: bkmr_", sex_level, "_", outcome_var, ".rda\n\n")
    }
  }
}
