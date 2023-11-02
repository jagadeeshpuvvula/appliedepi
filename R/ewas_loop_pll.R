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
  # Set the number of workers
  future::plan(future::multisession, workers = 20)

  #get lambda
  lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)

  # Iterate over exposures and bVals datasets
  for (variable in exposures) {
    for (i in seq_along(bVals_list)) {
      bVals <- bVals_list[[i]]
      bVals_name <- bVals_names[[i]]  # Get the dataset name

      if (variable %in% exposures) {
        message("Processing variable: ", variable, " with bVals dataset: ", bVals_name)

        ewas <- cpg.assoc(
          bVals,
          indep = dat_pheno[[variable]],
          covariates = as.data.frame(dat_pheno[, covariates]),
          logit.transform = TRUE)

        data_manhat <- cbind(ewas$results, ewas$coefficients) |> as_tibble()

        # Calculate and print lambda value
        lambda_val <- lambda(ewas$results[,3])
        message("Lambda value for variable ", variable, " with bVals dataset ", bVals_name, ": ", lambda_val)

        # Modify the file name for RDS to include bVals_name
        rds_file_name <- paste("ewas_", variable, "_", bVals_name, ".RData", sep = "")
        saveRDS(ewas, file = file.path(outputFolder, rds_file_name))

        # Modify the file name for CSV to include bVals_name
        csv_file_name <- paste("df_manhat_", variable, "_", bVals_name, ".csv", sep = "")
        write_csv(data_manhat, file = file.path(outputFolder, csv_file_name))

        message("Processing completed for variable: ", variable, " & methylation dataset: ", bVals_name)
      }
    }
  }
}
