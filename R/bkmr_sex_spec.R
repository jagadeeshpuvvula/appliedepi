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
                     iter=iter) {

  set.seed(2023)

  # Subset data based on sex level and update covariate list accordingly
  for (sex_level in c("male", "female")) {

    if (include_sex) {
      if (sex_level == "male") {
        sex_data <- subset(data, sex == "1")
        covariate_vars <- covariate_vars[!(covariate_vars %in% c("sex"))]
      } else if (sex_level == "female") {
        sex_data <- subset(data, sex == "2")
        covariate_vars <- covariate_vars[!(covariate_vars %in% c("sex"))]
      } else {
        stop("Invalid sex level")
      }
    } else {
      sex_data <- data
    }

    # Extract variables needed for bkmr
    ln_mixture <- as.matrix(sex_data[, ln_mixture_vars])
    outcome <- as.matrix(sex_data[, outcome_vars])
    covariates <- as.matrix(sex_data[, covariate_vars])

    # Run bkmr for each outcome variable
    for(i in 1:ncol(outcome)) {

      # Create knots for bkmr
      knots100 <- fields::cover.design(ln_mixture, nd = round(0.15 * nrow(ln_mixture)))$design

      bkmr <- kmbayes(y = outcome[,i], Z = ln_mixture, X = covariates, iter = iter,
                      verbose = TRUE, varsel = TRUE,
                      groups = groups,
                      knots = knots100)

      save(bkmr, file = paste0(folder_path, "/bkmr_",
                               if (include_sex) sex_level else "all", "_",
                               "all_",
                               colnames(outcome)[i], ".rda"))
    }

  }

}
