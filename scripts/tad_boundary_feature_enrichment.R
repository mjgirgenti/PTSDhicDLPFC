# ==============================================================================
# TAD boundary feature enrichment analysis
#
# This script tests whether differential TAD boundaries are enriched for
# genomic features relative to randomly sampled stable TAD boundaries.
#
# Feature enrichment is quantified as the proportion of boundaries overlapping
# at least one feature. For each differential boundary set, an empirical null
# distribution is generated from 1,000 random samples of stable boundaries
# containing the same number of boundaries as the observed set.
#
# Features tested:
#   - Gained DARs
#   - Lost DARs
#   - Strengthened chromatin loops
#   - Weakened chromatin loops
#   - Differentially expressed genes (DEGs)
#   - PTSD GWAS variants
#
# Required input objects:
#   - gained_dar_gr: GRanges of gained DARs
#   - lost_dar_gr: GRanges of lost DARs
#   - loop_anchor1_gr: GRanges of first anchors of differential loops
#   - loop_anchor2_gr: GRanges of second anchors of differential loops
#   - deg_gr: GRanges of DEGs
#   - gwas_gr: GRanges of PTSD GWAS variants
#   - altered_boundary_gr: GRanges of weakened or strengthened TAD boundaries
#   - stable_boundary_gr: GRanges of stable TAD boundaries
#
# Statistical analysis:
#   - Number of permutations: 1,000
#   - Random boundary sets are sampled without replacement from stable
#     boundaries and contain the same number of boundaries as the target set.
#   - Enrichment statistic: proportion of boundaries overlapping >=1 feature.
#   - Empirical enrichment P value: fraction of permutations with a proportion
#     >= the observed proportion.
#   - Empirical depletion P value: fraction of permutations with a proportion
#     <= the observed proportion.
# ==============================================================================


# ------------------------------------------------------------------------------
# R libraries
# ------------------------------------------------------------------------------

library(dplyr)
library(GenomicRanges)
library(pbapply)


# ==============================================================================
# Functions
# ==============================================================================


# ------------------------------------------------------------------------------
# Count unique chromatin loops overlapping each boundary
# ------------------------------------------------------------------------------

count_loops_per_boundary <- function(
    boundary_gr,
    anchor1_gr,
    anchor2_gr,
    direction = NULL
) {

  if (!is.null(direction)) {
    anchor1_use <- anchor1_gr[anchor1_gr$Direction == direction]
    anchor2_use <- anchor2_gr[anchor2_gr$Direction == direction]
  } else {
    anchor1_use <- anchor1_gr
    anchor2_use <- anchor2_gr
  }

  ov1 <- findOverlaps(
    boundary_gr,
    anchor1_use,
    ignore.strand = TRUE
  )

  ov2 <- findOverlaps(
    boundary_gr,
    anchor2_use,
    ignore.strand = TRUE
  )

  hit1 <- data.frame(
    boundary_id = queryHits(ov1),
    loop_id = anchor1_use$loop_id[subjectHits(ov1)]
  )

  hit2 <- data.frame(
    boundary_id = queryHits(ov2),
    loop_id = anchor2_use$loop_id[subjectHits(ov2)]
  )

  # Count each loop only once per boundary if both anchors overlap
  hits <- unique(rbind(hit1, hit2))

  counts <- integer(length(boundary_gr))

  if (nrow(hits) > 0) {
    tab <- table(hits$boundary_id)
    counts[as.integer(names(tab))] <- as.integer(tab)
  }

  return(counts)
}


# ------------------------------------------------------------------------------
# Count genomic features overlapping each boundary
# ------------------------------------------------------------------------------

count_features_per_boundary <- function(
    boundary_gr,
    gained_dar_gr,
    lost_dar_gr,
    loop_anchor1_gr,
    loop_anchor2_gr,
    deg_gr,
    gwas_gr,
    window_size = 40000
) {

  # Differential TAD boundaries are represented as 40-kb intervals.
  # For larger windows, resize each interval around its center.
  if (window_size == 40000) {
    boundary_window_gr <- boundary_gr
  } else {
    boundary_window_gr <- resize(
      boundary_gr,
      width = window_size,
      fix = "center"
    )
  }

  feature_counts <- data.frame(
    boundary_id = seq_along(boundary_window_gr),

    Gained_DAR = countOverlaps(
      boundary_window_gr,
      gained_dar_gr,
      ignore.strand = TRUE
    ),

    Lost_DAR = countOverlaps(
      boundary_window_gr,
      lost_dar_gr,
      ignore.strand = TRUE
    ),

    Strengthened_Loop = count_loops_per_boundary(
      boundary_window_gr,
      loop_anchor1_gr,
      loop_anchor2_gr,
      direction = "Strengthened"
    ),

    Weakened_Loop = count_loops_per_boundary(
      boundary_window_gr,
      loop_anchor1_gr,
      loop_anchor2_gr,
      direction = "Weakened"
    ),

    DEG = countOverlaps(
      boundary_window_gr,
      deg_gr,
      ignore.strand = TRUE
    ),

    GWAS = countOverlaps(
      boundary_window_gr,
      gwas_gr,
      ignore.strand = TRUE
    )
  )

  return(feature_counts)
}


# ------------------------------------------------------------------------------
# Randomly sample stable boundaries
# ------------------------------------------------------------------------------

sample_stable_boundaries <- function(
    stable_boundary_gr,
    n_boundaries
) {

  if (length(stable_boundary_gr) < n_boundaries) {
    stop(
      paste0(
        "Need ", n_boundaries,
        " stable boundaries, but only ",
        length(stable_boundary_gr), " are available."
      )
    )
  }

  sampled_idx <- sample(
    seq_along(stable_boundary_gr),
    size = n_boundaries,
    replace = FALSE
  )

  return(stable_boundary_gr[sampled_idx])
}


# ------------------------------------------------------------------------------
# Test feature enrichment at differential TAD boundaries
# ------------------------------------------------------------------------------

boundary_feature_enrichment <- function(
    altered_boundary_gr,
    stable_boundary_gr,
    gained_dar_gr,
    lost_dar_gr,
    loop_anchor1_gr,
    loop_anchor2_gr,
    deg_gr,
    gwas_gr,
    boundary_type,
    window_size = 40000,
    n_perm = 1000
) {

  feature_names <- c(
    "Gained_DAR",
    "Lost_DAR",
    "Strengthened_Loop",
    "Weakened_Loop",
    "DEG",
    "GWAS"
  )


  # --------------------------------------------------------------------------
  # Observed boundary-feature overlap
  # --------------------------------------------------------------------------

  observed_df <- count_features_per_boundary(
    boundary_gr = altered_boundary_gr,
    gained_dar_gr = gained_dar_gr,
    lost_dar_gr = lost_dar_gr,
    loop_anchor1_gr = loop_anchor1_gr,
    loop_anchor2_gr = loop_anchor2_gr,
    deg_gr = deg_gr,
    gwas_gr = gwas_gr,
    window_size = window_size
  )

  # Convert feature counts to presence/absence.
  # The statistic is the proportion of boundaries overlapping >=1 feature.
  observed_binary <- observed_df

  observed_binary[, feature_names] <-
    observed_binary[, feature_names, drop = FALSE] > 0

  observed_proportion <- colMeans(
    observed_binary[, feature_names, drop = FALSE]
  )


  # --------------------------------------------------------------------------
  # Empirical null distribution from stable boundaries
  # --------------------------------------------------------------------------

  random_matrix <- pbreplicate(
    n_perm,
    {

      sampled_stable_boundary_gr <- sample_stable_boundaries(
        stable_boundary_gr = stable_boundary_gr,
        n_boundaries = length(altered_boundary_gr)
      )

      random_df <- count_features_per_boundary(
        boundary_gr = sampled_stable_boundary_gr,
        gained_dar_gr = gained_dar_gr,
        lost_dar_gr = lost_dar_gr,
        loop_anchor1_gr = loop_anchor1_gr,
        loop_anchor2_gr = loop_anchor2_gr,
        deg_gr = deg_gr,
        gwas_gr = gwas_gr,
        window_size = window_size
      )

      random_df[, feature_names] <-
        random_df[, feature_names, drop = FALSE] > 0

      colMeans(
        random_df[, feature_names, drop = FALSE]
      )
    }
  )


  # --------------------------------------------------------------------------
  # Calculate fold enrichment and empirical P values
  # --------------------------------------------------------------------------

  result_df <- lapply(
    seq_along(feature_names),
    function(i) {

      feature <- feature_names[i]

      observed <- observed_proportion[[feature]]
      random_values <- random_matrix[i, ]

      random_count_high <- sum(random_values >= observed)
      random_count_low <- sum(random_values <= observed)

      # Empirical permutation P values used in the manuscript
      p_enrichment <- random_count_high / n_perm
      p_depletion <- random_count_low / n_perm

      data.frame(
        Boundary_type = boundary_type,
        Window_kb = window_size / 1000,
        Feature = feature,

        Observed_boundaries = sum(
          observed_binary[[feature]],
          na.rm = TRUE
        ),

        Observed_proportion = observed,

        Stable_mean_proportion = mean(random_values),
        Stable_SD = sd(random_values),

        Fold_enrichment = observed / mean(random_values),

        Random_count_high = random_count_high,
        Random_count_low = random_count_low,

        P_enrichment = p_enrichment,
        P_depletion = p_depletion,

        N_permutation = n_perm
      )
    }
  )

  return(bind_rows(result_df))
}



# ==============================================================================
# Example: NeuN+ weakened TAD boundaries, 40-kb window
# ==============================================================================

# ------------------------------------------------------------------------------
# Define input objects
# ------------------------------------------------------------------------------

# NP or NN
# Differential accessibility
gained_dar_gr <- gained_dar_gr
lost_dar_gr   <- lost_dar_gr

# Differential chromatin loop anchors
loop_anchor1_gr <- GRanges(
  seqnames = diff_loop_gr$chr,
  ranges = IRanges(
    start = diff_loop_gr$start1,
    end   = diff_loop_gr$end1))

loop_anchor2_gr <- GRanges(
  seqnames = diff_loop_gr$chr,
  ranges = IRanges(
    start = diff_loop_gr$start2,
    end   = diff_loop_gr$end2))
# add loop ID and Direction
loop_anchor1_gr$loop_id   <- seq_len(nrow(diff_loop_gr))
loop_anchor2_gr$loop_id   <- seq_len(nrow(diff_loop_gr))
loop_anchor1_gr$Direction <- diff_loop_gr$Direction
loop_anchor2_gr$Direction <- diff_loop_gr$Direction
loop_anchor1_gr <- loop_anchor1_gr
loop_anchor2_gr <- loop_anchor2_gr

# Differentially expressed genes
deg_gr <- deg_gr

# PTSD GWAS variants
gwas_gr <- PTSD_merged_snp38_sig_ld_uniq_gr

# TAD boundaries
altered_boundary_gr <- weakened_boundary_gr
stable_boundary_gr <- stable_boundary_gr




### test
# gained_dar_gr <- gained_dar_gr
# lost_dar_gr   <- lost_dar_gr
# loop_anchor1_gr <- loop_anchor1_gr
# loop_anchor2_gr <- loop_anchor2_gr
# deg_gr <- NP_deg_gr
# gwas_gr <- PTSD_merged_snp38_sig_ld_uniq_gr
# altered_boundary_gr <- NP_diff_tad_gr_weakened
# stable_boundary_gr <- res_coord_not_sig_gr



# ------------------------------------------------------------------------------
# Check input sizes
# ------------------------------------------------------------------------------

cat("NeuN+ 40-kb boundary enrichment input summary:\n")
cat("  Gained DARs:         ", length(gained_dar_gr), "\n")
cat("  Lost DARs:           ", length(lost_dar_gr), "\n")
cat("  Differential loops:  ", length(unique(loop_anchor1_gr$loop_id)), "\n")
cat("  DEGs:                ", length(deg_gr), "\n")
cat("  PTSD GWAS variants:  ", length(gwas_gr), "\n")
cat("  Target boundaries:   ", length(altered_boundary_gr), "\n")
cat("  Stable boundaries:   ", length(stable_boundary_gr), "\n")



# ------------------------------------------------------------------------------
# Run enrichment analysis
# ------------------------------------------------------------------------------

set.seed(123)

NP_weakened_40kb <- boundary_feature_enrichment(
  altered_boundary_gr = altered_boundary_gr,
  stable_boundary_gr = stable_boundary_gr,
  gained_dar_gr = gained_dar_gr,
  lost_dar_gr = lost_dar_gr,
  loop_anchor1_gr = loop_anchor1_gr,
  loop_anchor2_gr = loop_anchor2_gr,
  deg_gr = deg_gr,
  gwas_gr = gwas_gr,
  boundary_type = "Weakened",
  window_size = 40000,
  n_perm = 1000
)

NP_weakened_40kb





