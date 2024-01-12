#' match sequence - validation step before getting into MWAS loop
#' match the exposure dataset with the metabolome sequence
#' metabolome each observation is a feature and each variable is a subject
#' Whereas the exposure/covariate dataset each observation is a subject
#'
#' @param df1 A dataframe to match
#' @param df2 second dataframe to match with the df1
#' @import tidyverse
#' @export
#'
#' @return print's sequence of identifiers that matched between datasets
#' @examples
#' \donttest{
#' match_sequence(df1= mom_dat, #dataset with exposure/outcome variables
#'                df2=mom_ft, #dataset with metabolome features
#'                variable_name="participant_id_mom")
#' }
#'
match_sequence <- function(df1, df2, variable_name) {
  # Get the values of the specified variable in df1
  values_df1 <- df1[[variable_name]]

  # Get the variable names in df2
  variable_names_df2 <- names(df2)

  # Match the values of df1 to the variable names of df2
  matched_variable_names <- variable_names_df2[match(values_df1, variable_names_df2)]

  return(matched_variable_names)
}

