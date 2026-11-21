# ==============================================================================
# Differential chromatin accessibility analysis
#
# This script identifies differentially accessible regions (DARs) between
# PTSD and control samples using DESeq2. ATAC-seq peak counts are modeled while
# adjusting for demographic and ATAC-seq quality covariates.
#
# DARs are defined at FDR <= 0.05 and classified as gained or lost
# accessibility according to the direction of the log2 fold change.
#
# Required input objects:
#   - NeuN_pos_atac_count:
#       Matrix of ATAC-seq peak counts (peaks x samples) for NeuN+ nuclei.
#   - NeuN_pos_coldata:
#       Sample metadata corresponding to columns of the peak count matrix.
#   - normFactors:
#       Peak-specific normalization factors generated in the upstream
#       normalization analysis.
#
# Covariates:
#   - Sex
#   - Race
#   - AgeDeath
#   - PMI
#   - FRiP
#   - TSS_enrichment
# ==============================================================================


# ------------------------------------------------------------------------------
# R libraries
# ------------------------------------------------------------------------------

library(DESeq2)


# ------------------------------------------------------------------------------
# Prepare DESeq2 model
# ------------------------------------------------------------------------------

NeuN_pos_coldata$Condition <- factor(
  NeuN_pos_coldata$Condition,
  levels = c("Control", "Treatment"))
NeuN_pos_coldata$AgeDeath = scale(NeuN_pos_coldata$AgeDeath, center=TRUE)
NeuN_pos_coldata$PMI = scale(NeuN_pos_coldata$PMI, center=TRUE)
NeuN_pos_coldata$FRiP = scale(NeuN_pos_coldata$FRiP, center=TRUE)

dds <- DESeqDataSetFromMatrix(
  countData = NeuN_pos_atac_count,
  colData = NeuN_pos_coldata,
  design = ~ Sex +
    Race +
    AgeDeath +
    PMI +
    FRiP +
    TSS_enrichment +
    Condition
)


# ------------------------------------------------------------------------------
# Apply normalization factors
# ------------------------------------------------------------------------------

norm_factors_centered <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds) <- norm_factors_centered

# ------------------------------------------------------------------------------
# Likelihood ratio test
# ------------------------------------------------------------------------------

dds <- DESeq(
  dds,
  fitType = "glmGamPoi",
  test = "LRT",
  reduced = ~ Sex +
    Race +
    AgeDeath +
    PMI +
    FRiP +
    TSS_enrichment
)

results_all <- results(dds)


# ------------------------------------------------------------------------------
# Identify differentially accessible regions
# ------------------------------------------------------------------------------

differential_accessible_regions <- as.data.frame(
  results_all[
    !is.na(results_all$padj) &
      results_all$padj <= 0.05,
  ]
)

differential_accessible_regions$Direction <- ifelse(
  differential_accessible_regions$log2FoldChange > 0,
  "Gained",
  "Lost"
)


# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

dim(differential_accessible_regions)
table(differential_accessible_regions$Direction)