# ==============================================================================
# Differential TAD boundary analysis
#
# Differential TAD boundary strength was assessed using insulation scores.
# Insulation scores were quantile normalized across samples and analyzed using
# limma while adjusting for demographic and Hi-C quality covariates.
#
# Differential boundaries were defined using nominal P < 0.05.
#
# Required input objects:
#   - insulation_score_df:
#       Data frame containing genomic coordinates (chr, start, end) followed
#       by sample-level insulation scores.
#   - sample_info_hic:
#       Sample metadata with rows corresponding to columns of the insulation
#       score matrix.
#
# Covariates:
#   - Condition
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

library(limma)


# ------------------------------------------------------------------------------
# Prepare insulation score matrix
# ------------------------------------------------------------------------------

# Remove genomic coordinate columns (chr, start, end)
insulation_matrix <- as.matrix(
  insulation_score_df[, -c(1:3)]
)

# Sample order in the insulation matrix must match the metadata
stopifnot(
  identical(
    colnames(insulation_matrix),
    rownames(sample_info_hic)
  )
)


# ------------------------------------------------------------------------------
# Prepare sample metadata and design matrix
# ------------------------------------------------------------------------------

sample_info_hic$Condition <- factor(
  sample_info_hic$Condition,
  levels = c("Control", "Treatment")
)

sample_info_hic$Sex <- factor(sample_info_hic$Sex)

design <- model.matrix(
  ~ Condition +
    Sex +
    Race +
    AgeDeath +
    PMI +
    pct_valid +
    pct_valid_rmdup +
    pct_cis +
    pct_cis_short,
  data = sample_info_hic
)


# ------------------------------------------------------------------------------
# Quantile normalization
# ------------------------------------------------------------------------------

insulation_matrix_qn <- preprocessCore::normalize.quantiles(
  insulation_matrix,
  keep.names = TRUE
)


# ------------------------------------------------------------------------------
# Differential insulation score analysis
# ------------------------------------------------------------------------------

fit <- lmFit(
  insulation_matrix_qn,
  design
)

fit <- eBayes(fit)

results_all <- topTable(
  fit,
  coef = "ConditionTreatment",
  number = Inf
)


# Nominal P < 0.05 was used to identify differential boundaries
results_nominal_sig <- results_all[
  results_all$P.Value < 0.05,
]


# ------------------------------------------------------------------------------
# Add genomic coordinates to differential analysis results
# ------------------------------------------------------------------------------

# All tested TAD boundaries
all_boundary_df <- cbind(
  insulation_score_df[
    as.numeric(rownames(results_all)),
    c("chrom", "start", "end")
  ],
  results_all
)


# Differential TAD boundaries (nominal P < 0.05)
differential_boundary_df <- cbind(
  insulation_score_df[
    as.numeric(rownames(results_nominal_sig)),
    c("chrom", "start", "end")
  ],
  results_nominal_sig
)


# ------------------------------------------------------------------------------
# Define stable TAD boundaries
# ------------------------------------------------------------------------------

# Boundaries not meeting the nominal P < 0.05 threshold are operationally
# defined as stable boundaries for downstream permutation analyses.
stable_boundary_df <- insulation_score_df[
  -as.numeric(rownames(results_nominal_sig)),
  c("chrom", "start", "end")
]

# Add "chr" prefix if chromosome names are stored without it
stable_boundary_df$chrom <- ifelse(
  grepl("^chr", stable_boundary_df$chrom),
  stable_boundary_df$chrom,
  paste0("chr", stable_boundary_df$chrom)
)

stable_boundary_gr <- GRanges(
  stable_boundary_df
)
