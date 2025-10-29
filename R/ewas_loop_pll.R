#' Performs EWAS using CpGassoc library. Allows flexibility to loop for multiple exposure variables and methylation data from multiple tissues.
#'
#' @param exposures list of exposure variables
#' @param covariates list of covariates to adjust for associations
#' @param bVals_list processed beta values set from illumina methylation
#' @param bVals_names name of the beta value datasets - will be used for naming the output files
#' @param dat_pheno dataframe that contain exposures, covariates and cell type counts
#' @param outputFolder path to save the output files
#'
#' @return this function will return an R data object as a result of the CpG wide association using cpg.assoc function
#' @export
#'
#'@import CpGassoc
#'@import future.apply
#'@import future
#'
#' @example
#' \donttest{
#' ewas_loop_pl(exposures= c("fluo2", "nap1", "nap2", "phen1", "phen4", "phen9", "phen23", "pyr1"),
#'              covariates= c("mom_age", "mom_race", "mom_edu", "pgtob", "bmi",
#'              "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),
#'              bVals_list = list(bVals, mVals),
#'              bVals_names = list("x", "y"),
#'              dat_pheno = dat_pheno,
#'              outputFolder = "~/Documents/methylation_pilot/testing")
#'
#' }
ewas_loop_pll <- function(exposures, covariates, bVals_list, bVals_names, dat_pheno, outputFolder) {
  # Lambda calculation function
  lambda <- function(p) {
    p_clean <- p[!is.na(p)]
    if (length(p_clean) == 0) return(NA)
    median(qchisq(p_clean, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, df = 1)
  }
  
  # Iterate over exposures and bVals datasets
  for (variable in exposures) {
    for (i in seq_along(bVals_list)) {
      bVals <- bVals_list[[i]]
      bVals_name <- bVals_names[[i]]
      
      message("Processing variable: ", variable, " with bVals dataset: ", bVals_name)
      
      tryCatch({
        # Run EWAS
        ewas <- cpg.assoc(
          bVals,
          indep = dat_pheno[[variable]],
          covariates = as.data.frame(dat_pheno[, covariates]),
          logit.transform = TRUE
        )
        
        # Prepare results
        data_manhat <- cbind(ewas$results, ewas$coefficients) |> as_tibble()
        
        # Calculate lambda value with error handling
        lambda_val <- tryCatch(
          lambda(ewas$results[, 3]),
          error = function(e) {
            warning("Lambda calculation failed: ", e$message)
            NA
          }
        )
        message("Lambda value for variable ", variable, " with bVals dataset ", bVals_name, ": ", lambda_val)
        
        # Save results
        rds_file_name <- paste0("ewas_", variable, "_", bVals_name, ".rds")
        saveRDS(ewas, file = file.path(outputFolder, rds_file_name))
        
        csv_file_name <- paste0("df_manhat_", variable, "_", bVals_name, ".csv")
        write_csv(data_manhat, file = file.path(outputFolder, csv_file_name))
        
        message("Processing completed for variable: ", variable, " & methylation dataset: ", bVals_name)
        
      }, error = function(e) {
        warning("Error processing variable ", variable, " with dataset ", bVals_name, ": ", e$message)
      })
    }
  }
}
