#' CpG annotate and save csv files
#'
#' @param file_names list of csv file names that need to be annotated
#' @param input_folder location of the EWAS results in csv format
#' @param output_folder path to save the output files
#'
#' @return this function will return an R data object as a result of the CpG wide association using cpg.assoc function
#' @export
#'
#' @import tidyverse
#'
#' @example
#' \donttest{
#'annotate_and_save(file_names = c("df_manhat_fluo2_fp.csv",
#'                                 "df_manhat_fluo2_mp.csv",
#'                                 "df_manhat_pyr1_fp.csv",
#'                                 "df_manhat_pyr1_mp.csv",
#'                                 "df_manhat_nap1_mp.csv",
#'                                 "df_manhat_nap1_fp.csv",
#'                                 "df_manhat_nap1_mom_cmbc.csv",
#'                                 "df_manhat_phen4_mp.csv",
#'                                 "df_manhat_m_cnp_mp.csv",
#'                                 "df_manhat_bpa_mp.csv",
#'                                 "df_manhat_e_pb_mp.csv"),
#'                  input_folder = "~/Documents/methylation_pilot/res_cell_cnt",
#'                  output_folder = "~/Documents/methylation_pilot/res_annotated_cell_cnt")
#' }

annotate_and_save <- function(file_names, input_folder, output_folder) {
  # Get a list of all files in the input folder and its subfolders
  all_files <- list.files(path = input_folder, pattern = NULL, full.names = TRUE, recursive = TRUE)

  # Filter for the desired file names
  files_to_process <- all_files[basename(all_files) %in% file_names]

  for (file_path in files_to_process) {
    # Read the CSV file and clean column names
    df_res <- read_csv(file_path) |>
      clean_names()

    # Perform left join with ann450k data frame
    df_plot <- left_join(df_res, ann450k, by = "cpg_labels")

    # Create the output file path
    file_name <- basename(file_path)
    output_file_path <- file.path(output_folder, paste0("anno_", file_name))

    # Save the resulting data frame to the output file
    write_csv(df_plot, output_file_path)

    print(paste("Processed and saved", file_name, "to", output_file_path))
  }
}
