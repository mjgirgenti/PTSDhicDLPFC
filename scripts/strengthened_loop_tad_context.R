# ============================================================
# Strengthened loop analysis by TAD context in NeuN+ cells
# ============================================================
# This analysis tests whether the magnitude of PTSD-associated
# loop strengthening differs between loops contained within a
# single TAD and loops extending beyond individual TADs.
#
# Strengthened loops are classified according to whether the
# full genomic span between their two anchors is contained
# within a single TAD. Loop log2 fold changes are compared
# between the two groups, followed by linear regression
# adjusting for genomic loop span.
# ============================================================

library(GenomicRanges)
library(dplyr)



#loop_anchor1_gr or loop_anchor2_gr
# GRanges object with 708 ranges and 19 metadata columns:
#         seqnames              ranges strand |  baseMean log2FoldChange     lfcSE      stat      pvalue      padj   Direction
#            <Rle>           <IRanges>  <Rle> | <numeric>      <numeric> <logical> <numeric>   <numeric> <numeric> <character>
#     [1]    chr12     8100000-8110000      * |  12.08825        1.41117      <NA>   10.9671 0.000927445 0.0453705          Strengthened
#     [2]    chr13   92200000-92210000      * |   9.16614        1.66603      <NA>   12.4925 0.000408596 0.0323721          Weakened
#TAD_gr
# GRanges object with 5549 ranges and 3 metadata columns:
#          seqnames            ranges strand | boundary_type has_weakened_boundary has_strengthened_boundary
#             <Rle>         <IRanges>  <Rle> |      <factor>             <logical>                 <logical>
#      [1]     chr1   3880000-4640000      * |        Stable                 FALSE                     FALSE
#      [2]     chr1   4720000-5400000      * |        Stable                 FALSE                     FALSE


# ============================================================
# 1. Select strengthened loops
# ============================================================

strengthened_anchor1_gr <- loop_anchor1_gr[
    loop_anchor1_gr$log2FoldChange > 0
]

strengthened_anchor2_gr <- loop_anchor2_gr[
    loop_anchor2_gr$log2FoldChange > 0
]

stopifnot(
    length(strengthened_anchor1_gr) ==
        length(strengthened_anchor2_gr)
)

length(strengthened_anchor1_gr)
# 706


# ============================================================
# 2. Construct genomic spans between loop anchors
# ============================================================

loop_span_gr <- GRanges(
    seqnames = seqnames(strengthened_anchor1_gr),
    ranges = IRanges(
        start = pmin(
            start(strengthened_anchor1_gr),
            start(strengthened_anchor2_gr)
        ),
        end = pmax(
            end(strengthened_anchor1_gr),
            end(strengthened_anchor2_gr)
        )
    ),
    strand = "*"
)

length(loop_span_gr)
# 706


# ============================================================
# 3. Classify loops by TAD context
# ============================================================

hits_within <- findOverlaps(
    loop_span_gr,
    TAD_gr,
    type = "within",
    ignore.strand = TRUE
)

within_tad_idx <- unique(queryHits(hits_within))

within_tad <- seq_along(loop_span_gr) %in% within_tad_idx

table(within_tad)

# FALSE  TRUE
#   425   281

mean(within_tad)
# 0.398


# ============================================================
# 4. Construct analysis dataframe
# ============================================================

df_strength <- data.frame(
    log2FC = strengthened_anchor1_gr$log2FoldChange,
    TAD_context = factor(
        within_tad,
        levels = c(FALSE, TRUE),
        labels = c("Beyond TAD", "Within TAD")
    ),
    span = width(loop_span_gr)
)


# ============================================================
# 5. Summarize loop strengthening and genomic span
# ============================================================

strength_summary <- df_strength %>%
    group_by(TAD_context) %>%
    summarise(
        N = n(),
        median_log2FC = median(log2FC),
        mean_log2FC = mean(log2FC),
        median_span_kb = median(span) / 1000,
        mean_span_kb = mean(span) / 1000,
        .groups = "drop"
    )

strength_summary

#   TAD_context     N median_log2FC mean_log2FC median_span_kb mean_span_kb
# 1 Beyond TAD    425          1.63        1.67           700.         944.
# 2 Within TAD    281          1.37        1.39           330.         395.


# ============================================================
# 6. Compare strengthened-loop log2FC between TAD contexts
# ============================================================

wilcox_result <- wilcox.test(
    log2FC ~ TAD_context,
    data = df_strength
)

wilcox_result


# ============================================================
# 7. Linear regression adjusting for loop span
# ============================================================

fit <- lm(
    log2FC ~ TAD_context + log10(span),
    data = df_strength
)

summary(fit)
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -1.52555    0.24647  -6.190 1.03e-09 ***
# TAD_contextWithin TAD -0.09492    0.03039  -3.123  0.00186 ** 
# log10(span)            0.54616    0.04209  12.977  < 2e-16 ***




# ============================================================
# Strengthened loops by TAD boundary context in NeuN+ cells
# ============================================================
# Strengthened loops contained within a single TAD are further
# classified according to TAD boundary context (stable, weakened,
# strengthened, or mixed).
#
# Loop strengthening is compared across TAD contexts. Loop-span
# distributions are also tested to assess whether differences in
# strengthening could reflect systematic differences in genomic span.
# ============================================================


# ============================================================
# 1. Assign within-TAD strengthened loops to TAD context
# ============================================================

ov_loop_tad <- findOverlaps(
    loop_span_gr,
    TAD_gr,
    type = "within",
    ignore.strand = TRUE
)

loop_idx <- queryHits(ov_loop_tad)
tad_idx <- subjectHits(ov_loop_tad)

loop_tad_df <- data.frame(
    loop_id = loop_idx,
    tad_id = tad_idx,
    TAD_context = TAD_gr$boundary_type[tad_idx],
    log2FC = strengthened_anchor1_gr$log2FoldChange[loop_idx],
    span = width(loop_span_gr)[loop_idx]
)

table(loop_tad_df$TAD_context)
prop.table(table(loop_tad_df$TAD_context))

# Stable        248   (88.3%)
# Weakened       19   (6.8%)
# Strengthened   13   (4.6%)
# Mixed           1   (0.4%)


# ============================================================
# 2. Exclude mixed-boundary TADs
# ============================================================

loop_tad_test <- subset(
    loop_tad_df,
    TAD_context != "Mixed"
)

table(loop_tad_test$TAD_context)


# ============================================================
# 3. Compare loop spans across TAD contexts
# ============================================================

span_test <- pairwise.wilcox.test(
    log10(loop_tad_test$span),
    loop_tad_test$TAD_context,
    p.adjust.method = "BH",
    exact = FALSE
)

span_test

# BH-adjusted P:
# Weakened vs Stable        = 0.36
# Strengthened vs Stable    = 0.66
# Strengthened vs Weakened  = 0.66


# Median loop span by TAD context
span_summary <- loop_tad_test %>%
    group_by(TAD_context) %>%
    summarise(
        N = n(),
        median_span_kb = median(span) / 1000,
        .groups = "drop"
    )

span_summary

# Stable        330 kb
# Weakened      310 kb
# Strengthened  330 kb


# ============================================================
# 4. Compare loop strengthening across TAD contexts
# ============================================================

strength_test <- pairwise.wilcox.test(
    loop_tad_test$log2FC,
    loop_tad_test$TAD_context,
    p.adjust.method = "BH",
    exact = FALSE
)

strength_test

# BH-adjusted P:
# Weakened vs Stable        = 0.04
# Strengthened vs Stable    = 0.72
# Strengthened vs Weakened  = 0.04





