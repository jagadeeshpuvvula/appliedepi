#' Performs Pathway Enrichment Analysis (PEA) using mummichog function available from metaboAnalystR library
#' This function also prints the p-value theshold on scree that is being used for the analysis
#'
#' @param input_folder Specify folder location that contain the .txt files ready for the analysis. they should contain m.z, r.t, p.value variables
#' @param output_folder Specify location to save the output files
#' @param p_val_cutoff Specify a custom p-value or leave blank if "pick_500" is set to TRUE
#' @param pick_500 Set this parameter to TRUE if you need the function to pick a p-value threshold based on 500th smallest p-value
#'
#' @return returns mummichog results in csv format at the output_folder location
#' @export
#'
#' @examples
#' \donttest{
#'runMummichog(input_folder= "~/Documents/phth_phe_MWAS/result/enrich_exp",
#'             # p_val_cutoff= 0.2,
#'             pick_500 = TRUE,
#'             output_folder= "~/Documents/phth_phe_MWAS/result/enrich_imp")
#'
#' }
#'
#'
#'
runMummichog <- function(input_folder, output_folder, p_val_cutoff = NULL, pick_500 = FALSE) {
  # create object for storing data
  mSet3 <- InitDataObjects("mass_all", "mummichog", FALSE)

  # set peak format
  mSet3 <- SetPeakFormat(mSet3, "rmp")
  mSet3 <- UpdateInstrumentParameters(mSet3, 5.0, "mixed", "no")

  # get a list of text files in the input folder
  input_files <- list.files(input_folder, pattern = "\\.txt$", full.names = TRUE)

  # process each input file
  for (input_file in input_files) {
    # Get p_val_cutoff based on pick_500 parameter
    if (pick_500) {
      mummi_results <- read.table(input_file, header = TRUE)  # Assuming text files are tabular with headers
      if (nrow(mummi_results) >= 500) {
        sorted_pvals <- sort(mummi_results$p.value)
        p_val_cutoff <- sorted_pvals[500]
        cat("New p_val_cutoff (500th smallest p.val) for file", input_file, ":", p_val_cutoff, "\n")
      } else {
        cat("Data has less than 500 observations for file", input_file, ". Using input p_val_cutoff.\n")
        if (is.null(p_val_cutoff)) {
          stop("No p_val_cutoff provided and data has less than 500 observations.")
        }
      }
    } else {
      if (is.null(p_val_cutoff)) {
        stop("No p_val_cutoff provided and pick_500 is FALSE.")
      }
    }

    # read the peak list data
    mSet3 <- Read.PeakListData(mSet3, input_file)
    mSet3 <- SanityCheckMummichogData(mSet3)

    # map selected adducts to current data
    mSet3 <- Setup.AdductData(mSet3, add.vec)
    mSet3 <- PerformAdductMapping(mSet3, add.mode="mixed")

    # perform mummichog algorithm using selected adducts, using version 2 of the mummichog algorithm
    mSet3 <- SetPeakEnrichMethod(mSet3, algOpt="mum", version="v2")
    mSet3 <- SetMummichogPval(mSet3, cutoff=p_val_cutoff) # pval

    # the next step takes three or four minutes to run
    mSet3 <- PerformPSEA(mSet3, "hsa_mfn", "current", 3, 10000)

    # store the results as a data frame
    mummi_results <- as.data.frame(mSet3$mummi.resmat)

    # generate output file path and name
    output_file <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(input_file)), ".csv"))

    # save results as a CSV file
    write.csv(mummi_results, file = output_file, row.names = TRUE)

    # print the output file path for each processed file
    cat("Processed file:", input_file, "Output file:", output_file, "\n")
  }
}
