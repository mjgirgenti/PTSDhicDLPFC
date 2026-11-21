# ==============================================================================
# Differential chromatin loop analysis
#
# This script identifies differential chromatin loops between PTSD and control
# samples using DESeq2. Loop contact counts are modeled while adjusting for
# demographic and Hi-C quality covariates.
#
# Differential loops are defined at FDR <= 0.05 and classified as strengthened
# or weakened according to the direction of the log2 fold change.
#
# Required input objects:
#   - hic_loop_count_df:
#       Matrix of filtered loop contact counts (loops x samples).
#   - sample_info_hic:
#       Sample metadata corresponding to columns of the loop count matrix.
#   - normFactors:
#       Loop-specific normalization factors generated in the upstream
#       normalization analysis.
#
# Covariates:
#   - Sex
#   - Race
#   - AgeDeath
#   - PMI
#   - pct_valid
#   - pct_valid_rmdup
#   - pct_cis
#   - pct_cis_short
# ==============================================================================


# ------------------------------------------------------------------------------
# R libraries
# ------------------------------------------------------------------------------

library(DESeq2)


# ------------------------------------------------------------------------------
# Differential loop model
# ------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = hic_loop_count_df,
  colData = sample_info_hic,
  design = ~ Sex +
    Race +
    AgeDeath +
    PMI +
    pct_valid +
    pct_valid_rmdup +
    pct_cis +
    pct_cis_short +
    Condition
)


# ------------------------------------------------------------------------------
# Apply normalization factors
# ------------------------------------------------------------------------------

# Normalization factors
norm_factors_centered <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds) <- norm_factors_centered


# ------------------------------------------------------------------------------
# Likelihood ratio test
# ------------------------------------------------------------------------------
dds$Condition <- relevel(dds$Condition, ref = "Control")

dds <- DESeq(
  dds,
  fitType = "glmGamPoi",
  test = "LRT",
  reduced = ~ Sex +
    Race +
    AgeDeath +
    PMI +
    pct_valid +
    pct_valid_rmdup +
    pct_cis +
    pct_cis_short
)

results_all <- results(dds)


# ------------------------------------------------------------------------------
# Identify differential chromatin loops
# ------------------------------------------------------------------------------

differential_loops <- as.data.frame(
  results_all[
    !is.na(results_all$padj) &
      results_all$padj <= 0.05,
  ]
)

differential_loops$Direction <- ifelse(
  differential_loops$log2FoldChange > 0,
  "Strengthened",
  "Weakened"
)

table(differential_loops$Direction)

