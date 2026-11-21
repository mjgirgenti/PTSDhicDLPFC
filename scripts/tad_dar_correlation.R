# ============================================================
# Correlation between TAD boundary weakening and accessibility loss
# ============================================================
# This analysis tests whether the magnitude of TAD boundary
# weakening is associated with the magnitude of accessibility
# loss at lost DARs overlapping weakened TAD boundaries.
#
# For each weakened TAD boundary, the magnitudes of accessibility
# loss at overlapping lost DARs are averaged. Spearman correlation
# is then calculated between boundary weakening and mean
# accessibility loss.
#
# Inputs:
#   weakened_tad_boundaries_gr : weakened TAD boundaries with logFC metadata
#   lost_dar_gr     : lost DARs with log2FoldChange metadata
# ============================================================

library(GenomicRanges)
library(dplyr)


head(weakened_tad_boundaries_gr)
# GRanges object with 6 ranges and 6 metadata columns:
#        seqnames              ranges strand |     logFC    AveExpr         t     P.Value adj.P.Val         B
#           <Rle>           <IRanges>  <Rle> | <numeric>  <numeric> <numeric>   <numeric> <numeric> <numeric>
#   7348     chrX   41240000-41280000      * | 0.1754345 -0.2872619   4.69791 0.000716391  0.892804 -0.331619
#   4369     chr3   62280000-62320000      * | 0.1066159 -0.9120247   4.25029 0.001474472  0.892804 -0.900337

head(lost_dar_gr)
# GRanges object with 6 ranges and 18 metadata columns:
#        seqnames              ranges strand |    Length          Geneid log2FoldChange      pvalue        padj   Direction
#           <Rle>           <IRanges>  <Rle> | <integer>     <character>      <numeric>   <numeric>   <numeric> <character>
#   3529     chr7   56373408-56373893      + |       486 Interval_243471      -1.726273 8.94722e-38 2.83681e-33        Lost
#   5619    chr12 131576729-131577398      + |       670  Interval_67085      -4.186801 1.02739e-26 1.62872e-22        Lost


# ============================================================
# 1. Identify lost DARs overlapping weakened TAD boundaries
# ============================================================

ov <- findOverlaps(
    weakened_tad_boundaries_gr,
    lost_dar_gr,
    ignore.strand = TRUE
)

boundary_dar_df <- data.frame(
    boundary_id = queryHits(ov),
    dar_id = subjectHits(ov),

    # Magnitude of TAD boundary weakening
    boundary_weakening = abs(
        weakened_tad_boundaries_gr$logFC[queryHits(ov)]
    ),

    # Magnitude of accessibility loss
    accessibility_loss = abs(
        lost_dar_gr$log2FoldChange[subjectHits(ov)]
    )
)

head(boundary_dar_df)
nrow(boundary_dar_df)


# ============================================================
# 2. Summarize accessibility loss for each weakened boundary
# ============================================================

boundary_dar_fc <- boundary_dar_df %>%
    group_by(boundary_id) %>%
    summarise(
        boundary_weakening = unique(boundary_weakening),
        N_DAR = n(),
        mean_accessibility_loss = mean(accessibility_loss, na.rm = TRUE),
        .groups = "drop"
    )

# ============================================================
# 3. Test the association between boundary weakening
#    and accessibility loss
# ============================================================

cor_result <- cor.test(
    boundary_dar_fc$boundary_weakening,
    boundary_dar_fc$mean_accessibility_loss,
    method = "spearman",
    exact = FALSE
)

cor_result

# NeuN+:
# Spearman's rho = -0.60
# P = 0.025
#
# NeuN-:
# Spearman's rho = -0.05
# P = 0.837


