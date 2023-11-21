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
#' \donttest{
#' process_methylation_data(basenames= idat_files)
#' }
process_methylation_data <- function(basenames) {
  # Read methylation data
  rgSet <- read.metharray(basenames = basenames)

  # Print initial dimensions
  cat("Initial dimensions: ", nrow(rgSet), " rows, ", ncol(rgSet), " columns\n")

  # Remove samples based on detection p-values
  #Drop samples if 5% of probes are above the p-value threshold
  detP <- detectionP(rgSet)
  samples_to_remove <- colSums(detP > 1e-7) > (0.05 * nrow(rgSet))
  rgset_f2 <- rgSet[!samples_to_remove, ]

  # Print dimensions after sample removal
  cat("Dimensions after sample removal based on detection p-values: ", nrow(rgset_f2), " rows, ", ncol(rgset_f2), " columns\n")
  cat("Total rows dropped: ", nrow(rgSet) - nrow(rgset_f2), "\n")
  cat("Total columns dropped: ", ncol(rgSet) - ncol(rgset_f2), "\n")

  # Preprocess with noob and bmiq
  mSetSq_noob <- preprocessNoob(rgset_f2)
  mSetSq_bmiq <- BMIQ(mSetSq_noob)

  # Batch correction with combat
  mSetSq_combat <- ComBat(dat = mSetSq_bmiq, batch = pheno_by_tiss$batch)

  # Convert preprocessed +BMIQ object to a genomic ratio set
  mSetSq_v1 <- makeGenomicRatioSetFromMatrix(mSetSq_combat)

  # Remove CpGs that have failed in one or more samples
  keep <- rowSums(detP < 1e-7) == ncol(mSetSq_v1)
  mSetSqFlt <- mSetSq_v1[keep,]

  # Print dimensions after CpG removal
  cat("Dimensions after removing CpGs that have failed in one or more samples: ", nrow(mSetSqFlt), " rows, ", ncol(mSetSqFlt), " columns\n")
  cat("Total rows dropped: ", nrow(mSetSq_v1) - nrow(mSetSqFlt), "\n")
  cat("Total columns dropped: ", ncol(mSetSq_v1) - ncol(mSetSqFlt), "\n")

  # Remove CpGs on the sex chromosomes
  keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX", "chrY")])
  mSetSqFlt <- mSetSqFlt[keep,]

  # Print dimensions after removing CpGs on the sex chromosomes
  cat("Dimensions after removing CpGs on the sex chromosomes: ", nrow(mSetSqFlt), " rows, ", ncol(mSetSqFlt), " columns\n")

  # Remove CpG's with cross-reactive probes
  xReactiveProbes <- read_csv(url("https://epigen.ccm.sickkids.ca/sample-report/data/quality_control/cross_reactive_probes.csv"),
                              col_names = c("TargetID"))
  keep <- !(rownames(mSetSqFlt) %in% xReactiveProbes$TargetID)
  mSetSqFlt <- mSetSqFlt[keep,]

  # Print dimensions after removing CpGs with cross-reactive probes
  cat("Dimensions after removing CpGs with cross-reactive probes: ", nrow(mSetSqFlt), " rows, ", ncol(mSetSqFlt), " columns\n")

  # Remove probes with SNPs at CpG site
  mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

  # Print dimensions after removing probes with SNPs at CpG site
  cat("Dimensions after removing probes with SNPs at CpG site: ", nrow(mSetSqFlt), " rows, ", ncol(mSetSqFlt), " columns\n")

  return(mSetSqFlt)
}
