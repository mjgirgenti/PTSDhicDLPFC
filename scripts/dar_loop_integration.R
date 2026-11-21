# ==============================================================================
# DAR–differential loop integration
#
# This analysis tests whether gained or lost DARs preferentially overlap
# anchors of strengthened or weakened chromatin loops.
#
# Required input objects:
#   - differential_loop_gi:
#       GRanges containing differential loop anchors. Each anchor contains
#       a "Direction" field indicating strengthened or weakened loops.
#
#   - dar_gr:
#       GRanges containing DARs. Each region contains a "Direction" field
#       indicating gained or lost accessibility.
#
# Statistical analysis:
#   - DARs overlapping at least one differential loop anchor are identified
#     using genomic overlap.
#   - Fisher's exact test compares the proportions of gained and lost DARs
#     overlapping loop anchors.
# ==============================================================================


library(GenomicRanges)


head(differential_loop_gi)
# GInteractions object with 6 interactions and 6 metadata columns:
#       seqnames1             ranges1     seqnames2             ranges2 |  baseMean log2FoldChange      stat      pvalue      padj    Direction
#           <Rle>           <IRanges>         <Rle>           <IRanges> | <numeric>      <numeric> <numeric>   <numeric> <numeric>  <character>
#   [1]     chr12     8100000-8110000 ---     chr12     8370000-8380000 |  12.08825        1.41117   10.9671 0.000927445 0.0453705 Strengthened
#   [2]     chr13   92200000-92210000 ---     chr13   92660000-92670000 |   9.16614        1.66603   12.4925 0.000408596 0.0323721 Strengthened

head(dar_gr)
# GRanges object with 6 ranges and 6 metadata columns:
#       seqnames            ranges strand |  baseMean log2FoldChange      stat      pvalue        padj   Direction
#          <Rle>         <IRanges>  <Rle> | <numeric>      <numeric> <numeric>   <numeric>   <numeric> <character>
#   [1]     chr8 12192993-12195168      * |  620.1309       1.592074   213.553 6.08785e-27 7.47624e-22      Gained
#   [2]    chr21   8390074-8391024      * | 3589.0180      -0.805315   145.542 1.78593e-21 7.31077e-17        Lost


# ------------------------------------------------------------------------------
# Function to count DARs overlapping differential loop anchors
# ------------------------------------------------------------------------------

count_loop_overlapping_dar <- function(
    loop_gr,
    dar_gr,
    loop_direction,
    dar_direction
) {

  loop_subset <- loop_gr[
    loop_gr$Direction == loop_direction
  ]

  dar_subset <- dar_gr[
    dar_gr$Direction == dar_direction
  ]

  # Identify DARs overlapping either anchor of differential chromatin loops
  overlaps <- findOverlaps(
    loop_subset,
    dar_subset,
    ignore.strand = TRUE
  )

  # Count each DAR only once even if it overlaps multiple loop anchors
  length(unique(subjectHits(overlaps)))
}


# ------------------------------------------------------------------------------
# Count overlaps between DARs and differential loop anchors
# ------------------------------------------------------------------------------

n_gained_strengthened <- count_loop_overlapping_dar(
  loop_gr = differential_loop_gi,
  dar_gr = dar_gr,
  loop_direction = "Strengthened",
  dar_direction = "Gained"
)

n_lost_strengthened <- count_loop_overlapping_dar(
  loop_gr = differential_loop_gi,
  dar_gr = dar_gr,
  loop_direction = "Strengthened",
  dar_direction = "Lost"
)

n_gained_weakened <- count_loop_overlapping_dar(
  loop_gr = differential_loop_gi,
  dar_gr = dar_gr,
  loop_direction = "Weakened",
  dar_direction = "Gained"
)

n_lost_weakened <- count_loop_overlapping_dar(
  loop_gr = differential_loop_gi,
  dar_gr = dar_gr,
  loop_direction = "Weakened",
  dar_direction = "Lost"
)

c(
  Gained_Strengthened = n_gained_strengthened,
  Lost_Strengthened = n_lost_strengthened,
  Gained_Weakened = n_gained_weakened,
  Lost_Weakened = n_lost_weakened
)


# Gained_Strengthened    37
# Lost_Strengthened      16
# Gained_Weakened         0
# Lost_Weakened           0


# ------------------------------------------------------------------------------
# Fisher's exact test: strengthened loops
# ------------------------------------------------------------------------------

n_gained_dar <- sum(
  dar_gr$Direction == "Gained"
)

n_lost_dar <- sum(
  dar_gr$Direction == "Lost"
)

loop_dar_table <- matrix(
  c(
    n_gained_strengthened,
    n_lost_strengthened,
    n_gained_dar - n_gained_strengthened,
    n_lost_dar - n_lost_strengthened
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Loop_overlap = c("Yes", "No"),
    DAR_direction = c("Gained", "Lost")
  )
)

loop_dar_table

loop_dar_fisher <- fisher.test(
  loop_dar_table
)

loop_dar_fisher
#p-value = 1.199e-09



