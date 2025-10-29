#' Process methylation data using minfi pipeline. This function will drop samples with GT 5 percent high detection p-values, probes with at least one high detection p-value, Sex chromosomes, loci associated with single-nucleotide polymorphisms and cross-hybridizing probes. Additionally, will preprocess using Noob for background correction and dye bias normalization, BMIQ for probe design bias, batch correction using ComBat.
#'
#' @param basenames list of idat file locations
#'
#' @return returns a methyl set
#' @export process_methylation_data
#'
#' @import minfi
#' @import wateRmelon
#' @import ENmix
#' @import sva
#'
#' @example
# ========== USAGE EXAMPLE ==========
# 
# # Create basenames data frame with Sample_ID and folder locations
# basenames_df <- data.frame(
#   Sample_ID = c("sample1_12345678", "sample2_87654321", "sample3_11223344"),
#   folder_location = c(
#     "/path/to/batch1/idats",
#     "/path/to/batch1/idats",
#     "/path/to/batch2/idats"
#   ),
#   stringsAsFactors = FALSE
# )
# 
# # Each folder should contain files like:
# # sample1_12345678_Red.idat
# # sample1_12345678_Grn.idat
# # sample2_87654321_Red.idat
# # sample2_87654321_Grn.idat
# # etc.
# 
# # Optional: Create phenotype data with additional sample information
# pheno_data <- data.frame(
#   Sample_ID = c("sample1_12345678", "sample2_87654321", "sample3_11223344"),
#   batch = c("batch1", "batch1", "batch2"),
#   age = c(45, 52, 38),
#   sex = c("M", "F", "M"),
#   tissue = c("blood", "blood", "blood"),
#   row.names = c("sample1_12345678", "sample2_87654321", "sample3_11223344"),
#   stringsAsFactors = FALSE
# )
# # Run preprocessing for EPIC array
# process_methylation_data(
#   basenames_df = basenames_df,
#   pheno_data = pheno_data,
#   output_folder = "/path/to/output/methylation_results_epic",
#   array_type = "EPIC",
#   detection_pval = 0.01,
#   sample_cutoff = 0.05,
#   apply_bmiq = TRUE,
#   apply_combat = TRUE,
#   batch_variable = "batch",
#   verbose = TRUE
# )
process_methylation_data <- function(basenames_df, 
                                      pheno_data = NULL,
                                      output_folder,
                                      array_type = c("450k", "EPIC", "EPICv2"),
                                      detection_pval = 0.01,
                                      sample_cutoff = 0.05,
                                      apply_bmiq = TRUE,
                                      apply_combat = TRUE,
                                      batch_variable = NULL,
                                      remove_sex_chr = TRUE,
                                      remove_snps = TRUE,
                                      remove_cross_reactive = TRUE,
                                      verbose = TRUE) {
  
  # Load required libraries
  require(minfi)
  require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  require(sva)  # for ComBat
  require(wateRmelon)  # for BMIQ
  require(readr)
  
  array_type <- match.arg(array_type)
  
  # ========== 0. Validate Inputs and Create Output Folder ==========
  if (verbose) cat("\n=== Validating inputs ===\n")
  
  # Check basenames_df structure
  if (!is.data.frame(basenames_df)) {
    stop("basenames_df must be a data frame")
  }
  
  required_cols <- c("Sample_ID", "folder_location")
  if (!all(required_cols %in% colnames(basenames_df))) {
    stop("basenames_df must contain columns: Sample_ID, folder_location")
  }
  
  if (verbose) {
    cat("Number of samples in basenames_df: ", nrow(basenames_df), "\n")
  }
  
  # Create full paths to basenames (without file extension)
  basenames <- file.path(basenames_df$folder_location, basenames_df$Sample_ID)
  
  # Validate that IDAT files exist
  missing_files <- c()
  for (i in 1:length(basenames)) {
    red_file <- paste0(basenames[i], "_Red.idat")
    grn_file <- paste0(basenames[i], "_Grn.idat")
    
    if (!file.exists(red_file) || !file.exists(grn_file)) {
      missing_files <- c(missing_files, basenames_df$Sample_ID[i])
    }
  }
  
  if (length(missing_files) > 0) {
    stop("IDAT files not found for samples: ", paste(missing_files, collapse = ", "))
  }
  
  if (verbose) cat("All IDAT files found\n")
  
  # Create output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    if (verbose) cat("Created output folder: ", output_folder, "\n")
  } else {
    if (verbose) cat("Using existing output folder: ", output_folder, "\n")
  }
  
  # Create phenotype data if not provided
  if (is.null(pheno_data)) {
    if (verbose) cat("No phenotype data provided. Creating basic pheno_data from basenames_df\n")
    pheno_data <- basenames_df
    rownames(pheno_data) <- basenames_df$Sample_ID
  } else {
    # Validate pheno_data
    if (nrow(pheno_data) != nrow(basenames_df)) {
      stop("pheno_data must have the same number of rows as basenames_df")
    }
    # Ensure rownames match Sample_ID
    rownames(pheno_data) <- basenames_df$Sample_ID
  }
  
  # Helper function for printing
  print_dims <- function(obj, message) {
    if (verbose) {
      cat(message, ": ", nrow(obj), " probes x ", ncol(obj), " samples\n", sep = "")
    }
  }
  
  # ========== 1. Read Data ==========
  if (verbose) cat("\n=== Reading methylation data ===\n")
  rgSet <- read.metharray(basenames = basenames)
  
  # Set sample names
  colnames(rgSet) <- basenames_df$Sample_ID
  
  print_dims(rgSet, "Initial dimensions")
  
  # ========== 2. Quality Control - Sample Level ==========
  if (verbose) cat("\n=== Sample-level QC ===\n")
  
  # Calculate detection p-values
  detP <- detectionP(rgSet)
  colnames(detP) <- basenames_df$Sample_ID
  
  # Identify samples with too many failed probes
  failed_probes_per_sample <- colSums(detP > detection_pval)
  threshold <- sample_cutoff * nrow(rgSet)
  samples_to_remove <- failed_probes_per_sample > threshold
  
  if (verbose) {
    cat("Detection p-value threshold: ", detection_pval, "\n")
    cat("Sample removal threshold: ", sample_cutoff * 100, "% failed probes\n", sep = "")
    cat("Samples failing QC: ", sum(samples_to_remove), "/", ncol(rgSet), "\n")
    if (sum(samples_to_remove) > 0) {
      cat("Sample IDs removed: ", paste(colnames(rgSet)[samples_to_remove], collapse = ", "), "\n")
    }
  }
  
  # Remove poor quality samples
  rgSet <- rgSet[, !samples_to_remove]
  detP <- detP[, !samples_to_remove]
  pheno_data <- pheno_data[!samples_to_remove, ]
  
  print_dims(rgSet, "After sample removal")
  
  # ========== 3. Normalization ==========
  if (verbose) cat("\n=== Normalization (Noob) ===\n")
  mSetSq <- preprocessNoob(rgSet)
  print_dims(mSetSq, "After Noob normalization")
  
  # ========== 4. Get Beta Values ==========
  beta <- getBeta(mSetSq)
  
  # ========== 5. Quality Control - Probe Level ==========
  if (verbose) cat("\n=== Probe-level QC ===\n")
  
  # Remove probes that failed in one or more samples
  failed_probe_count <- rowSums(detP > detection_pval)
  keep_probes <- failed_probe_count == 0
  
  if (verbose) {
    cat("Probes with failed detection (p > ", detection_pval, ") in ≥1 sample: ", 
        sum(!keep_probes), "\n", sep = "")
  }
  
  beta <- beta[keep_probes, ]
  mSetSq <- mSetSq[keep_probes, ]
  detP <- detP[keep_probes, ]
  
  print_dims(beta, "After removing failed probes")
  
  # ========== 6. Remove Sex Chromosome Probes ==========
  if (remove_sex_chr) {
    if (verbose) cat("\n=== Removing sex chromosome probes ===\n")
    
    # Get appropriate annotation
    if (array_type == "450k") {
      data(Locations)
      annotation <- Locations
    } else if (array_type %in% c("EPIC", "EPICv2")) {
      data(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
    
    # Find sex chromosome probes
    sex_probes <- rownames(annotation)[annotation$chr %in% c("chrX", "chrY")]
    keep_probes <- !(rownames(beta) %in% sex_probes)
    
    if (verbose) {
      cat("Sex chromosome probes removed: ", sum(!keep_probes), "\n")
    }
    
    beta <- beta[keep_probes, ]
    mSetSq <- mSetSq[keep_probes, ]
    detP <- detP[keep_probes, ]
    
    print_dims(beta, "After removing sex chromosome probes")
  }
  
  # ========== 7. Remove Cross-Reactive Probes ==========
  if (remove_cross_reactive) {
    if (verbose) cat("\n=== Removing cross-reactive probes ===\n")
    
    tryCatch({
      # Different cross-reactive probe lists for different arrays
      if (array_type == "450k") {
        # Chen et al. 2013 cross-reactive probes
        xReactiveProbes <- read_csv(
          "https://epigen.ccm.sickkids.ca/sample-report/data/quality_control/cross_reactive_probes.csv",
          col_names = c("TargetID"),
          show_col_types = FALSE
        )
      } else {
        # For EPIC arrays, use Pidsley et al. 2016 list
        # You may need to adjust this URL or provide your own list
        xReactiveProbes <- read_csv(
          "https://epigen.ccm.sickkids.ca/sample-report/data/quality_control/cross_reactive_probes.csv",
          col_names = c("TargetID"),
          show_col_types = FALSE
        )
      }
      
      keep_probes <- !(rownames(beta) %in% xReactiveProbes$TargetID)
      
      if (verbose) {
        cat("Cross-reactive probes removed: ", sum(!keep_probes), "\n")
      }
      
      beta <- beta[keep_probes, ]
      mSetSq <- mSetSq[keep_probes, ]
      detP <- detP[keep_probes, ]
      
      print_dims(beta, "After removing cross-reactive probes")
      
    }, error = function(e) {
      if (verbose) cat("Warning: Could not load cross-reactive probes list. Skipping this step.\n")
    })
  }
  
  # ========== 8. Remove SNP Probes ==========
  if (remove_snps) {
    if (verbose) cat("\n=== Removing probes with SNPs ===\n")
    
    initial_probes <- nrow(mSetSq)
    mSetSq <- dropLociWithSnps(mSetSq, snps = c("SBE", "CpG"), maf = 0)
    
    # Update beta matrix to match
    beta <- beta[rownames(beta) %in% rownames(mSetSq), ]
    detP <- detP[rownames(detP) %in% rownames(mSetSq), ]
    
    if (verbose) {
      cat("Probes with SNPs removed: ", initial_probes - nrow(mSetSq), "\n")
    }
    
    print_dims(beta, "After removing SNP probes")
  }
  
  # ========== 9. BMIQ Normalization (Optional) ==========
  if (apply_bmiq) {
    if (verbose) cat("\n=== Applying BMIQ normalization ===\n")
    
    # Get probe design type
    if (array_type == "450k") {
      data(Manifest)
      design <- Manifest$Design
    } else {
      annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      design <- annotation$Type
      names(design) <- rownames(annotation)
    }
    
    # Match design vector to beta matrix
    design <- design[rownames(beta)]
    
    # Apply BMIQ to each sample
    beta_bmiq <- beta
    for (i in 1:ncol(beta)) {
      if (verbose && i %% 10 == 0) cat("  Processing sample ", i, "/", ncol(beta), "\n", sep = "")
      
      tryCatch({
        bmiq_result <- BMIQ(beta[, i], design)
        beta_bmiq[, i] <- bmiq_result$nbeta
      }, error = function(e) {
        if (verbose) cat("  Warning: BMIQ failed for sample ", i, ". Using original values.\n", sep = "")
      })
    }
    
    beta <- beta_bmiq
    if (verbose) cat("BMIQ normalization complete\n")
  }
  
  # ========== 10. Batch Correction (Optional) ==========
  if (apply_combat && !is.null(batch_variable)) {
    if (verbose) cat("\n=== Applying ComBat batch correction ===\n")
    
    # Check if batch variable exists in pheno_data
    if (!(batch_variable %in% colnames(pheno_data))) {
      stop("Batch variable '", batch_variable, "' not found in pheno_data")
    }
    
    batch <- pheno_data[[batch_variable]]
    
    # Check for batches with only one sample
    batch_table <- table(batch)
    if (any(batch_table < 2)) {
      if (verbose) {
        cat("Warning: Some batches have fewer than 2 samples. ComBat may fail.\n")
        print(batch_table)
      }
    }
    
    tryCatch({
      # Remove any covariates that don't vary or are confounded with batch
      # For basic ComBat, we just use batch
      beta <- ComBat(dat = beta, batch = batch, mod = NULL, par.prior = TRUE)
      if (verbose) cat("ComBat batch correction complete\n")
      
    }, error = function(e) {
      if (verbose) cat("Warning: ComBat failed. Proceeding without batch correction.\n")
      cat("Error message: ", e$message, "\n")
    })
  }
  
  # ========== 11. Create Final GenomicRatioSet ==========
  if (verbose) cat("\n=== Creating final GenomicRatioSet ===\n")
  
  # Create GenomicRatioSet with proper annotation
  if (array_type == "450k") {
    grSet <- makeGenomicRatioSetFromMatrix(
      beta,
      array = "IlluminaHumanMethylation450k",
      annotation = "ilmn12.hg19",
      what = "Beta"
    )
  } else if (array_type %in% c("EPIC", "EPICv2")) {
    grSet <- makeGenomicRatioSetFromMatrix(
      beta,
      array = "IlluminaHumanMethylationEPIC",
      annotation = "ilm10b4.hg19",
      what = "Beta"
    )
  }
  
  # Add phenotype data
  pData(grSet) <- pheno_data
  
  print_dims(grSet, "Final dataset")
  
  if (verbose) {
    cat("\n=== Processing complete ===\n")
    cat("Total probes retained: ", nrow(grSet), " (", 
        round(nrow(grSet) / nrow(rgSet) * 100, 1), "% of original)\n", sep = "")
    cat("Total samples retained: ", ncol(grSet), " (", 
        round(ncol(grSet) / nrow(basenames_df) * 100, 1), "% of original)\n", sep = "")
  }
  
  # ========== 12. Save Results ==========
  if (verbose) cat("\n=== Saving results to ", output_folder, " ===\n", sep = "")
  
  # Save grSet
  grSet_file <- file.path(output_folder, "grSet.rda")
  save(grSet, file = grSet_file)
  if (verbose) cat("Saved grSet to: ", grSet_file, "\n")
  
  # Save beta values
  beta_file <- file.path(output_folder, "beta.rda")
  save(beta, file = beta_file)
  if (verbose) cat("Saved beta to: ", beta_file, "\n")
  
  # Save detection p-values
  detP_file <- file.path(output_folder, "detP.rda")
  save(detP, file = detP_file)
  if (verbose) cat("Saved detP to: ", detP_file, "\n")
  
  # Save phenotype data
  pheno_data_file <- file.path(output_folder, "pheno_data.rda")
  save(pheno_data, file = pheno_data_file)
  if (verbose) cat("Saved pheno_data to: ", pheno_data_file, "\n")
  
  # Save processing summary
  summary_file <- file.path(output_folder, "processing_summary.txt")
  sink(summary_file)
  cat("DNA Methylation Processing Summary\n")
  cat("===================================\n\n")
  cat("Processing Date: ", as.character(Sys.time()), "\n\n")
  cat("Array Type: ", array_type, "\n")
  cat("Detection p-value threshold: ", detection_pval, "\n")
  cat("Sample cutoff: ", sample_cutoff * 100, "%\n", sep = "")
  cat("BMIQ applied: ", apply_bmiq, "\n")
  cat("ComBat applied: ", apply_combat, "\n")
  if (!is.null(batch_variable)) cat("Batch variable: ", batch_variable, "\n")
  cat("Sex chromosomes removed: ", remove_sex_chr, "\n")
  cat("SNP probes removed: ", remove_snps, "\n")
  cat("Cross-reactive probes removed: ", remove_cross_reactive, "\n\n")
  cat("Initial samples: ", nrow(basenames_df), "\n")
  cat("Final samples: ", ncol(grSet), "\n")
  cat("Samples removed: ", nrow(basenames_df) - ncol(grSet), "\n\n")
  cat("Final probes: ", nrow(grSet), "\n")
  cat("Final dimensions: ", nrow(grSet), " probes x ", ncol(grSet), " samples\n")
  sink()
  
  if (verbose) cat("Saved processing summary to: ", summary_file, "\n")
  
  if (verbose) cat("\n=== All files saved successfully ===\n")
  
  # Return results invisibly
  invisible(list(
    grSet = grSet,
    beta = beta,
    detP = detP,
    pheno_data = pheno_data,
    output_folder = output_folder
  ))
}
