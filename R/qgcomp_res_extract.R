#' Extract results from the qgcomp output RDA files. Iterates though the RDA files and compiles results in tibble format
#'
#' @param results_folder input a folder location that contain qgcomp RDA output files
#'
#' @return results a dataframe that will contain bootstap or not, gender, cohort, outcome, estimate, confidence intervals
#' @export
#'
#' @example
#' \donttest{
#'results<- get_gcomp_estimates(results_folder = "E:/BBK17/pj/data_2023apr/results/qgcomp/qgcomp_pcs_wo_wisc")
#' }
#'
#'
#'
#'
get_gcomp_estimates <- function(results_folder) {
    # get the file names in the folder
    file_names <- list.files(results_folder, pattern = "\\.rda$", full.names = FALSE)
    
    # create empty data frame to store the results
    results <- data.frame(boot_strp = character(),
                          gender = character(),
                          cohort = character(),
                          outcome = character(),
                          estimate = numeric(),
                          lower_ci = numeric(),
                          upper_ci = numeric(),
                          p_value = numeric(),
                          sample_size = numeric(),
                          stringsAsFactors = FALSE)
    
    # loops through all .rda files and extract estimates, CI & p-values
    for (file_name in file_names) {
        # load the object from the file
        file_path <- file.path(results_folder, file_name)
        
        # Get objects before loading
        objs_before <- ls()
        load(file_path)
        objs_after <- ls()
        
        # Find the newly loaded object
        new_objs <- setdiff(objs_after, objs_before)
        object_name <- new_objs[grep("^boot|^nb", new_objs)]
        
        if (length(object_name) == 0) {
            warning(paste("No object matching pattern found in file:", file_name))
            next
        }
        
        if (length(object_name) > 1) {
            warning(paste("Multiple objects found in file:", file_name, "- using first one"))
            object_name <- object_name[1]
        }
        
        # Get the object
        obj <- get(object_name)
        
        # extract the required information
        estimate <- obj$coef["psi1"]
        CI <- obj$ci
        lower_ci <- CI[1]
        upper_ci <- CI[2]
        p_value <- obj$pval[2]
        
        # Extract sample size from qx dimensions
        sample_size <- if (!is.null(obj$qx) && !is.null(dim(obj$qx))) {
            dim(obj$qx)[1]
        } else if (!is.null(obj$qx) && is.vector(obj$qx)) {
            length(obj$qx)
        } else {
            # Debug: print structure to see what's available
            cat("  Available components:", names(obj), "\n")
            NA
        }
        
        # split file name by "_" and remove ".rda"
        file_parts <- strsplit(gsub("\\.rda$", "", file_name), "_")[[1]]
        
        # Check if we have enough parts
        if (length(file_parts) < 5) {
            warning(paste("File name does not have expected format:", file_name))
            next
        }
        
        # assign parts to variables
        part1 <- file_parts[1]
        part2 <- file_parts[2]
        part3 <- file_parts[3]
        part4 <- file_parts[5]
        
        # print out the values
        cat("File:", file_name, "\n")
        cat("  Estimate:", estimate, "\n")
        cat("  Lower CI:", lower_ci, "\n")
        cat("  Upper CI:", upper_ci, "\n")
        cat("  p-value:", p_value, "\n")
        cat("  Sample size:", sample_size, "\n")
        cat("  Part 1:", part1, "\n")
        cat("  Part 2:", part2, "\n")
        cat("  Part 3:", part3, "\n")
        cat("  Part 4:", part4, "\n")
        
        # add the results to the data frame
        results <- rbind(results, data.frame(boot_strp = part1,
                                             gender = part2,
                                             cohort = part3,
                                             outcome = part4,
                                             estimate = estimate,
                                             lower_ci = lower_ci,
                                             upper_ci = upper_ci,
                                             p_value = p_value,
                                             sample_size = sample_size,
                                             stringsAsFactors = FALSE))
        
        # remove the object from the R environment to avoid conflicts
        rm(list = object_name)
    }
    
    return(results)
}
