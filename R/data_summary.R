#' Summarize categorical and continuous variables
#'
#' @param data input a dataframe that contain variables
#' @param continuous_variables provide list of categorical variables
#' @param categorical_variables provide list of continuous variables
#'
#' @return this function will return an object with 2 lists. one contain summary of continuous and other with categorical variables
#' @export
#'
#' @example
#' \donttest{
#'summaries_overall <- summarize_variables(
#'                          continuous_variables= c("cotinine", "mom_age", "home_score_total",
#'                                                "gest_age", "b_wght", "b_length"),
#'                          categorical_variables= c("sex", "parity_n", "edu3", "fish_int_pm", "race"),
#'                          data= dat)
#' }
summarize_variables <- function(data, continuous_variables, categorical_variables) {
  # Convert categorical variables to factors if they are not already
  data[categorical_variables] <- lapply(data[categorical_variables], as.factor)

  # Summary for categorical variables
  summarize_categorical <- function(data, categorical_vars) {
    variable_names <- c()
    summary_values <- c()

    for (var in categorical_vars) {
      counts <- table(data[[var]], useNA = "always")
      percentages <- prop.table(counts) * 100
      percentages_rounded <- sprintf("%.2f", round(percentages, 2))

      summary_string <- paste(counts, "(",percentages_rounded,")")

      variable_names <- c(variable_names, paste(var, names(counts), sep = "_"))
      summary_values <- c(summary_values, summary_string)
    }

    categorical_summary_df <- data.frame(Variable = variable_names, Value = summary_values)
    return(categorical_summary_df)
  }

  # Summary for continuous variables
  summarize_continuous <- function(data, continuous_vars) {
    variable_names <- c()
    summary_values <- c()

    for (var in continuous_vars) {
      na_count <- sum(is.na(data[[var]]))
      na_percent <- sprintf("%.2f", sum(is.na(data[[var]])) / length(data[[var]]) * 100)
      median_val <- sprintf("%.2f", round(median(data[[var]], na.rm = TRUE), 2))
      q1_val <- sprintf("%.2f", round(quantile(data[[var]], 0.25, na.rm = TRUE), 2))
      q3_val <- sprintf("%.2f", round(quantile(data[[var]], 0.75, na.rm = TRUE), 2))

      summary_string <- paste("NA:", na_count, "(",na_percent,")",
                              "Median:", median_val,
                              "(",q1_val, ", ", q3_val,")"
      )

      variable_names <- c(variable_names, var)
      summary_values <- c(summary_values, summary_string)
    }

    continuous_summary_df <- data.frame(Variable = variable_names, Value = summary_values)
    return(continuous_summary_df)
  }

  # Generating summaries
  categorical_summary <- summarize_categorical(data, categorical_variables)
  continuous_summary <- summarize_continuous(data, continuous_variables)

  # Return summaries in separate data frames
  return(list(categorical_summary = categorical_summary, continuous_summary = continuous_summary))
}
