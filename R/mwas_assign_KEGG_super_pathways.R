#' Assigns KEGG super pathways to the mummichog results (load curated list first)
#'
#' @param data input dataframe that contain all mummichog results
#' @param group_column dummy variable that takes the group names
#' @param pathway_column variable that has the pathway names
#' @param pathway_list load this from the separate R file available in this library
#'
#' @return
#' @export
#'
#' @examples
assignGroup <- function(data, group_column, pathway_column, pathway_list) {
  group_column <- "other"

  for (i in 1:nrow(data)) {
    pathway <- tolower(data[[pathway_column]][i])

    min_distance <- Inf
    closest_match <- ""

    for (group_name in names(pathway_list)) {
      group_data <- lapply(pathway_list[[group_name]], tolower)

      # Calculate string distances and find the closest match
      distances <- sapply(group_data, function(x) adist(x, pathway, partial = TRUE))
      min_dist <- min(distances)

      if (min_dist < min_distance) {
        min_distance <- min_dist
        closest_match <- group_name
      }
    }

    # If the closest match is below a certain threshold, assign the group
    if (min_distance < 3) {  # Adjust the threshold as needed
      data$group_column[i] <- closest_match
    }
  }

  return(data)
}
