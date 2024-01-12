#'Manually get p-value cutoffs based on the 500 features with lowest p-values
#'use along with the runMummichog function
#'
#'
#' @param input_file
#' @param default_cutoff
calculatePvalCutoff <- function(input_file, default_cutoff) {
  mummi_results <- read.table(input_file, header = TRUE)  # Assuming text files are tabular with headers
  if (nrow(mummi_results) >= 500) {
    sorted_pvals <- sort(mummi_results$p.value)
    p_val_cutoff <- sorted_pvals[500]
    cat("New p_val_cutoff (500th smallest p.val) for file", input_file, ":", p_val_cutoff, "\n")
    return(p_val_cutoff)
  } else {
    cat("Data has less than 500 observations for file", input_file, ". Using default p_val_cutoff.\n")
    return(default_cutoff)
  }
}
