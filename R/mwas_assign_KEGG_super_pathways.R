#' Add KEGG superpathways to the pathways enriched from the PEA analysis
#' Run pathway_list before using this function or will result in errors
#'
#' @param data dataframe that contain mummichog results in global environment
#' @param pathway_column name of the variable that has pathway names
#' @param pathway_list Curated list of KEGG pathways for mapping pathways to groups
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' dat <- assignGroup(dat,
#'                     pathway_column = "pathway",
#'                     pathway_list = pathways)
#' }
#'
#'
#'
assignGroup <- function(data, pathway_column, pathway_list) {
  data$group <- "other"

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
      data$group[i] <- closest_match
    }
  }

  return(data)
}
