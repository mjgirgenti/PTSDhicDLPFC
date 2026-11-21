# ============================================================
# TAD boundary hotspot analysis
# ============================================================
# The genome is divided into 1-Mb bins, and weakened TAD
# boundaries overlapping each bin are counted separately for
# NeuN+ and NeuN- cells.
#
# Genome-wide background boundary density (lambda) is defined
# as the total number of weakened boundaries divided by the
# total number of 1-Mb genomic bins.
#
# Boundary clustering is evaluated under a Poisson model using
# the upper-tail probability P(X >= k). A stringent significance
# threshold of P < 0.001 is used to define boundary enrichment.
# In both NeuN+ and NeuN-, >=3 weakened boundaries per 1-Mb bin
# satisfies this threshold and is therefore used to define
# weakened TAD boundary hotspots.
# ============================================================

head(weakened_tad_boundaries_gr)
# GRanges object with 6 ranges and 6 metadata columns:
#        seqnames              ranges strand |     logFC    AveExpr         t     P.Value adj.P.Val         B
#           <Rle>           <IRanges>  <Rle> | <numeric>  <numeric> <numeric>   <numeric> <numeric> <numeric>
#   7348     chrX   41240000-41280000      * | 0.1754345 -0.2872619   4.69791 0.000716391  0.892804 -0.331619
#   4369     chr3   62280000-62320000      * | 0.1066159 -0.9120247   4.25029 0.001474472  0.892804 -0.900337


library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# ============================================================
# 1. Generate 1-Mb genomic bins
# ============================================================

window_size <- 1e6
keep_chr <- paste0("chr", c(1:22, "X", "Y"))

windows <- tileGenome(
    seqinfo(txdb),
    tilewidth = window_size,
    cut.last.tile.in.chrom = TRUE
)

windows <- windows[seqnames(windows) %in% keep_chr]


# ============================================================
# 2. Count weakened TAD boundaries in each 1-Mb bin
# ============================================================

windows$boundary_count <- countOverlaps(
    windows,
    weakened_tad_boundaries_gr
)

table(windows$boundary_count)

windows[
    order(windows$boundary_count, decreasing = TRUE)
][1:10]


# ============================================================
# 3. Calculate genome-wide background boundary density
# ============================================================

total_boundaries <- length(weakened_tad_boundaries_gr)
total_bins <- length(windows)

lambda <- total_boundaries / total_bins

lambda


# ============================================================
# 4. Evaluate Poisson probabilities for boundary clustering
# ============================================================

poisson_tail <- function(k, lambda) {
    ppois(
        k - 1,
        lambda = lambda,
        lower.tail = FALSE
    )
}

max_count <- max(windows$boundary_count)

poisson_results <- data.frame(
    boundary_count = 2:max_count,
    p_value = sapply(
        2:max_count,
        poisson_tail,
        lambda = lambda
    )
)

poisson_results

# ============================================================
# 5. Define weakened TAD boundary hotspots
# ============================================================

#NeuN+
hotspot_threshold <- 3
#NeuN-
#hotspot_threshold <- 4

windows$hotspot <- windows$boundary_count >= hotspot_threshold

table(windows$hotspot)
sum(windows$hotspot)


# ============================================================
# 6. Extract hotspot bins
# ============================================================

hotspot_gr <- windows[windows$hotspot]

hotspot_gr <- hotspot_gr[
    order(hotspot_gr$boundary_count, decreasing = TRUE)
]

hotspot_gr


