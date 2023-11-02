#' Gene enrichment for EWAS files
#'
#' @param input_folder Annotated version of EWAS results in csv format
#' @param output_folder location to save the results
#'
#' @return returns a csv that will contain pathways associated with input CpGs with q-value <0.05 as threshold
#' @export
#'
#' @example
#' \donttest{
perform_gene_enrichment(input_folder = "~/Documents/methylation_pilot/res_annotated_cell_cnt",
                        output_folder = "~/Documents/methylation_pilot/res_enriched_cell_cnt")
#' }
perform_gene_enrichment <- function(input_folder, output_folder) {
  # Get a list of files in the input folder
  input_files <- list.files(path = input_folder, full.names = TRUE)

  # Check if the output folder exists, if not, create it
  if (!file.exists(output_folder)) {
    dir.create(output_folder)
  }

  # Iterate over each input file
  for (input_file in input_files) {
    # Load the input file (assuming it's in CSV format, modify as needed)
    df <- read.csv(input_file)

    # Perform gene enrichment
    gene_enri <- gometh(
      sig.cpg = df$cpg_labels[df$fdr < 0.05],
      all.cpg = df$cpg_labels,
      collection = "KEGG"
    ) |> rownames_to_column("pathway_id")

    # Extract the filename without the prefix "anno_df_manhat_"
    filename <- gsub("^anno_df_manhat_", "", basename(input_file))

    # Save the output to the output folder without the prefix
    output_file <- file.path(output_folder, paste0(filename))
    write_csv(gene_enri, file = output_file)

    # Print a message indicating completion for this file
    cat("Gene enrichment completed for", basename(input_file), "\n")
  }
}
