#' Prints detailed summary of the R environment including dependencies
#'
#' @return displays R environment summary
#' @export
#'
#' @examples
#' \donttest{
#' print_version_info()
#' }
print_version_info <- function() {
  # Print R version information
  cat( R.version$version.string, "\n\n")

  # Get the names of loaded libraries
  loaded_libs <- search()[grepl("^package:", search())]
  loaded_libs <- gsub("^package:", "", loaded_libs)

  # Print loaded library versions
  cat("Dependencies:\n")
  for (lib in loaded_libs) {
    version <- packageVersion(lib)
    cat(paste(lib, ":", version, "\n"))
  }
}
