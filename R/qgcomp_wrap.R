#' Quantile g-computation wrapper function to automate analysis for multiple outcomes. Additionally, allows flexibility to perform stratified analysis.
#'
#' @param outcomes input a list of outcome variables
#' @param dat input dataset in tibble format
#' @param output_folder complete location to a folder to save the output RDA files
#' @param chemicals list of variables that need to be considered as exposure mixture
#' @param covariates list of covariate variables
#' @param include_sex input TRUE or FALSE for sex specific stratification
#' @param include_cohort input TRUE or FALSE for cohort specific stratification
#' @param q input a numeric value for number of quantiles
#' @param b input a numeric value for number of random bootstrap samples
#' @param seed input a numeric value for random sample generation for bootstrap process
#'
#' @return This function will return RDA output files for each model run by boot and non-bootstrap method. Filenames will specify boot or nonboot, sex or cohort specific, and outcome name.
#' @export
#'
#'@import qgcomp
#'
#' @example
#' \donttest{
#'qgcomp_func(outcomes = c("wppsi_fsiq", "wppsi_viq", "wppsi_piq"),
#'            data = dat,
#'            output_folder = "E:/BBK17/pj/data_2023apr/results/qgcomp/qgcomp_29chem",
#'            include_sex = TRUE,
#'            include_cohort = TRUE,
#'            chemicals = list("log_Pb", "log_Hg", "log_DMA", "log_DDE", "log_PBDE_47", "log_PCB_118", "log_PCB_138", "log_PCB_153", "log_PCB_180", "log_PFHxS", "log_PFOA", "log_PFOS", "log_BCEtP", "log_BDCPP", "log_DBuP", "log_DPhP", "log_TCS", "log_BPA", "log_MBP", "log_MBZP", "log_MCPP", "log_sigma_DEHP", "log_MEP", "log_MIBP", "log_di_Ethyl_OP", "log_di_Methyl_OP", "log_B_PB", "log_M_PB", "log_P_PB"),
#'            covariates = list("race_bin", "log_cotinine", "mom_edu_cat", "home_score_total", "parity_n", "mom_age", "fish_int_pm"),
#'            q = 10, b=400, seed = 2022)
#'
#' }
#'
qgcomp_func <- function(outcomes, data, output_folder, chemicals, covariates,
                        include_sex = TRUE, include_cohort = TRUE, q, b, seed) {

  sex_levels <- if (include_sex) c("all", "Female", "Male") else "all"
  cohort_levels <- if (include_cohort) c("all", "home", "mirec") else "all"

  for (sex_level in sex_levels) {
    if (sex_level == "all") {
      sex_data <- data
      sex_formula <- "sex"
    } else {
      sex_data <- subset(data, sex == sex_level)
      sex_formula <- ""
    }

    for (cohort_level in cohort_levels) {
      if (cohort_level == "all") {
        cohort_data <- sex_data
        cohort_formula <- "+cohort + city"
        filename_prefix <- paste0(sex_level, "_all")
      } else if (cohort_level == "home") {
        cohort_data <- subset(sex_data, cohort == "1")
        cohort_formula <- ""
        filename_prefix <- paste0(sex_level, "_home")
      } else if (cohort_level == "mirec") {
        cohort_data <- subset(sex_data, cohort == "2")
        cohort_formula <- "+city"
        filename_prefix <- paste0(sex_level, "_mirec")
      } else {
        stop("Invalid cohort level")
      }

      for (outcome in outcomes) {
        formula <- as.formula(paste(outcome, "~",
                                    paste(chemicals, collapse = " + "), "+",
                                    sex_formula, "+",
                                    cohort_formula, "+",
                                    paste(covariates, collapse = " + ")))

        nb <- qgcomp.noboot(formula, expnms = mixture, data = cohort_data, family = gaussian(), q = q)
        boot <- qgcomp.boot(formula, expnms = mixture, data = cohort_data, family = gaussian(), q = q, B = b, seed = seed)
        save(nb, file = paste0(output_folder, "/", "nb_", filename_prefix, "_", outcome, ".rda"))
        save(boot, file = paste0(output_folder, "/", "boot_", filename_prefix, "_", outcome, ".rda"))
      }
    }
  }
}
