% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mwas_lm_loop.R
\name{mwas_lm_loop}
\alias{mwas_lm_loop}
\title{mwas loop - perform linear regression between each feature and exposure metric provided in the input
metabolome data will be transposed in the function
input metabolome features and exposures are in different dataframes
runtime 17,000 * 70 * 17 matrix takes around 10 mins}
\usage{
mwas_lm_loop(
  feature_table,
  exp_cov_data,
  output_folder,
  mwas_file_name,
  cutoff_file_name,
  exposures,
  covar,
  fdr_cutoff
)
}
\arguments{
\item{feature_table}{provide a metabolome dataframe where each observation is a metabolome feature and each column is a subject}

\item{exp_cov_data}{each row represent a subject and contains exposure and covariate variables as columns}

\item{output_folder}{provide the location to the folder to save the output file from this function}

\item{mwas_file_name}{name of the CSV file with the linear regression results that will be save in output_folder}

\item{exposures}{provide the variable names that represent the exposures}

\item{covar}{provide list of covariates that need to be adjusted in the model}

\item{fdr_cutoff}{specify the FDR cutoff value here}
}
\value{
This function will return two CSV files in the output_folder location provided by user 1st file contain linear regression results and 2nd with FDR threshold values
}
\description{
mwas loop - perform linear regression between each feature and exposure metric provided in the input
metabolome data will be transposed in the function
input metabolome features and exposures are in different dataframes
runtime 17,000 * 70 * 17 matrix takes around 10 mins
}
\examples{
\donttest{
mwas_lm_loop(feature_table = mom_ft,
            exp_cov_data = mom_dat,
            exposures = names(mom_dat)[4:20],
            covar = c("bmi", "mom_age", "pgtob", "mom_edu", "mom_race"),
            output_folder = "~/Documents/phth_phe_MWAS/result/lm_mwas",
            mwas_file_name = "mom_mwas.csv", #compiles all results in a single csv
            fdr_cutoff = 0.2)
}
}
