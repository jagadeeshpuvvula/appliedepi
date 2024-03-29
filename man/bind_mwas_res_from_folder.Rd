% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mwas_bind_csv_files.R
\name{bind_mwas_res_from_folder}
\alias{bind_mwas_res_from_folder}
\title{This function will bind rows across the csv files in a given target folder}
\usage{
bind_mwas_res_from_folder(folder_path)
}
\arguments{
\item{folder_path}{User need to prove the folder location that contain csv files that end with mwas
these csv files with mwas as a prefix will be resulted from the mwas_lm_loop function above
Especially helpful if the results are in multiple files}
}
\value{
this function will return a dataframe within the global environment
}
\description{
This function will bind rows across the csv files in a given target folder
}
\examples{
\donttest{
bind_mwas_res_from_folder(paste0(result,"lm_mwas"))
}


}
