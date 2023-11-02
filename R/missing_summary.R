#' Missing data summary, prints frequency and count
#'
#' @param df input a dataframe
#'
#' @import tidyverse
#' @import scales
#'
#' @return returns a dataframe that contain frequency and percent of missing observations
#' @export
#'
#' @example
#' \donttest{
#' dat_miss<- missing_data_summary(dat_basc_all_exp)
#' }
#'
#'
missing_data_summary <- function(df) {
  # Create a dataframe with the count of missing values for each variable
  missing_counts <- df |>
    summarise_all(~ sum(is.na(.))) |>
    gather(variable, missing_count)

  # Filter out variables with no missing values
  missing_counts <- missing_counts |>
    filter(missing_count > 0)

  # Calculate the percent of missing values for each variable
  missing_percents <- df |>
    summarise_all(~ mean(is.na(.))) |>
    gather(variable, missing_percent)

  # Merge the count and percent dataframes
  missing_summary <- left_join(missing_counts, missing_percents, by = "variable")

  # Print the summary table
  missing_summary |>
    mutate(missing_percent = scales::percent(missing_percent)) |>
    arrange(desc(missing_count))
}
