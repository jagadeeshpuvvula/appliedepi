#' Multiple linear regression wrapper function that will take a list of exposure and outcome variables. Returns a matrix of results between each combination
#'
#' @param dependent_vars input a list of outcome variables
#' @param independent_vars input a list of exposure variables
#' @param covariates input a list of confounders
#' @param data input a dataset tested with dataframe and tibble formats
#' @param include_sex input either TRUE or FALSE to perform stratified analysis
#' @param include_cohort input either TRUE or FALSE to perform stratified analysis
#' @param conf_level input confidence interval level. Default to yield 95 percent Confidence interval
#'
#' @return This function will return a dataframe that will contain exposure, outcome, cohort level, sex level, regression coefficient, 95 percent confidence interval
#' @export
#'
#' @example
#'\donttest{
#'lm_func(dependent_vars = c("wppsi_fsiq", "wppsi_viq", "wppsi_piq"),
#'                  independent_vars=c("log_Pb", "log_Hg", "log_DMA", "log_DDE", "log_PBDE_47"),
#'                  covariates = c("log_cotinine", "home_score_total", "mom_age"),
#'                  dat = dat, include_sex = TRUE, include_cohort = TRUE)
#'}
#'
#'@import stats
#'
lm_func <- function(dependent_vars, independent_vars, data, covariates, include_sex = TRUE, include_cohort = TRUE, conf_level = 0.95) {
  # Create empty lists to store results
  dependent_list <- list()
  independent_list <- list()
  sex_level_list <- list()
  cohort_level_list <- list()
  coef_list <- list()
  p_value_list <- list()
  ci_lower_list <- list()
  ci_upper_list <- list()
  n_obs_list <- list()  # NEW: for observation count
  
  if ("sex" %in% colnames(data)) {
    sex_present <- TRUE
  } else {
    sex_present <- FALSE
    include_sex <- FALSE
  }
  
  if ("cohort" %in% colnames(data)) {
    cohort_present <- TRUE
  } else {
    cohort_present <- FALSE
    include_cohort <- FALSE
  }
  
  sex_levels <- if (include_sex) c("all", "Female", "Male") else "all"
  cohort_levels <- if (include_cohort) c("all", "home", "mirec") else "all"
  
  for (sex_level in sex_levels) {
    if (sex_level == "all" || !sex_present) {
      sex_data <- data
      sex_formula <- ""
    } else {
      sex_data <- subset(data, sex == sex_level)
      sex_data <- droplevels(sex_data)
      sex_formula <- ""
    }
    
    for (cohort_level in cohort_levels) {
      if (cohort_level == "all" || !cohort_present) {
        cohort_data <- sex_data
        cohort_formula <- ""
      } else if (cohort_level == "home") {
        cohort_data <- subset(sex_data, cohort == "1")
        cohort_data <- droplevels(cohort_data)
        cohort_formula <- ""
      } else if (cohort_level == "mirec") {
        cohort_data <- subset(sex_data, cohort == "2")
        cohort_data <- droplevels(cohort_data)
        cohort_formula <- "+city"
      } else {
        stop("Invalid cohort level")
      }
      
      # Loop through all combinations of dependent and independent variables
      for (i in 1:length(dependent_vars)) {
        for (j in 1:length(independent_vars)) {
          # Create the covariates formula dynamically based on the covariates argument
          covariate_formula <- paste(covariates, collapse = " + ")
          
          # Run linear regression with dynamic covariates
          formula <- as.formula(paste(dependent_vars[i], "~", independent_vars[j], "+", covariate_formula, sex_formula, cohort_formula))
          
          # Fit model with na.action = na.exclude to handle NAs
          model <- lm(formula, data = cohort_data, na.action = na.omit)
          
          # Calculate confidence interval
          ci <- confint(model, parm = independent_vars[j], level = conf_level)
          
          # Store results in lists
          dependent_list[[length(dependent_list) + 1]] <- dependent_vars[i]
          independent_list[[length(independent_list) + 1]] <- independent_vars[j]
          sex_level_list[[length(sex_level_list) + 1]] <- sex_level
          cohort_level_list[[length(cohort_level_list) + 1]] <- cohort_level
          coef_list[[length(coef_list) + 1]] <- coef(model)[independent_vars[j]]
          p_value_list[[length(p_value_list) + 1]] <- summary(model)$coefficients[independent_vars[j], 4]
          ci_lower_list[[length(ci_lower_list) + 1]] <- ci[1]
          ci_upper_list[[length(ci_upper_list) + 1]] <- ci[2]
          n_obs_list[[length(n_obs_list) + 1]] <- nobs(model)  # NEW: number of observations
        }
      }
    }
  }
  
  # Create a dataframe with results
  results <- data.frame(
    dependent_variable = unlist(dependent_list),
    independent_variable = unlist(independent_list),
    sex_level = unlist(sex_level_list),
    cohort_level = unlist(cohort_level_list),
    coefficient = unlist(coef_list),
    p_value = unlist(p_value_list),
    ci_lower = unlist(ci_lower_list),
    ci_upper = unlist(ci_upper_list),
    n_obs = unlist(n_obs_list)  # NEW: add observation count
  )
  
  # Add FDR-adjusted p-values (Benjamini-Hochberg method)
  # Calculate FDR within each sex_level and cohort_level combination
  results <- results |>
    group_by(sex_level, cohort_level) |>
    mutate(p_value_fdr = p.adjust(p_value, method = "fdr")) |>
    ungroup()
  
  # Return the dataframe
  return(results)
}
