# ============================================================
# DAR target-gene and GRE analysis
# ============================================================
# This analysis identifies genes associated with PTSD-related
# DARs through promoter overlap and downstream regulatory
# annotations, and evaluates glucocorticoid response element
# (GRE)-containing DARs.
#
# DARs are first annotated relative to gene features using
# ChIPseeker. Gained and lost DARs are analyzed separately.
# Promoter-associated genes are defined from DARs annotated
# within +/-3 kb of transcription start sites.
#
# Input:
#   dar_gr : significant DARs with Direction metadata
# ============================================================

library(ChIPseeker)
library(GenomicRanges)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


head(dar_gr)
# GRanges object with 6 ranges and 6 metadata columns:
#       seqnames            ranges strand |  baseMean log2FoldChange      stat      pvalue        padj   Direction
#          <Rle>         <IRanges>  <Rle> | <numeric>      <numeric> <numeric>   <numeric>   <numeric> <character>
#   [1]     chr8 12192993-12195168      * |  620.1309       1.592074   213.553 6.08785e-27 7.47624e-22      Gained
#   [2]    chr21   8390074-8391024      * | 3589.0180      -0.805315   145.542 1.78593e-21 7.31077e-17        Lost


# ============================================================
# 1. Annotate DARs relative to gene features
# ============================================================

dar_anno <- annotatePeak(
    dar_gr,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
)

dar_anno_gr <- dar_anno@anno

length(dar_anno_gr)


# ============================================================
# 2. Separate gained and lost DARs
# ============================================================

gained_dar_gr <- dar_anno_gr[
    dar_anno_gr$Direction == "Gained"
]

lost_dar_gr <- dar_anno_gr[
    dar_anno_gr$Direction == "Lost"
]

length(gained_dar_gr)
length(lost_dar_gr)


# ============================================================
# 3. Identify promoter-associated genes
# ============================================================

gained_promoter_genes <- unique(
    gained_dar_gr$SYMBOL[
        grepl("Promoter", gained_dar_gr$annotation)
    ]
)

gained_promoter_genes <- na.omit(gained_promoter_genes)

length(gained_promoter_genes)
#[1] 807

lost_promoter_genes <- unique(
    lost_dar_gr$SYMBOL[
        grepl("Promoter", lost_dar_gr$annotation)
    ]
)

lost_promoter_genes <- na.omit(lost_promoter_genes)

length(lost_promoter_genes)
#[1] 806




# ============================================================
# Loop-based target-gene assignment for DARs
# ============================================================
# This analysis identifies promoter-associated target genes of
# DARs through chromatin loops.
#
# A loop-linked target gene is defined when:
#   1. one loop anchor overlaps a DAR;
#   2. the opposite anchor does not overlap a DAR; and
#   3. the opposite anchor overlaps a promoter (+/-3 kb from TSS).
#
# Both loop orientations are considered:
#   DAR at anchor 1 -> promoter at anchor 2
#   DAR at anchor 2 -> promoter at anchor 1
#
# Inputs:
#   dar_gr       : DARs used for target-gene assignment
#   loop_anchor1 : first loop anchors
#   loop_anchor2 : second loop anchors
# ============================================================

library(ChIPseeker)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


head(gained_dar_gr)
# GRanges object with 6 ranges and 0 metadata columns:
#        seqnames            ranges strand
#           <Rle>         <IRanges>  <Rle>
#   1246    chr21   8390074-8391024      +
#   2870     chr8 57209517-57210369      +

head(loop_anchor1)
# GRanges object with 6 ranges and 0 metadata columns:
#         seqnames              ranges strand
#            <Rle>           <IRanges>  <Rle>
#     [1]    chr12     8100000-8110000      *
#     [2]    chr13   92200000-92210000      *
head(loop_anchor2)
# GRanges object with 6 ranges and 0 metadata columns:
#       seqnames              ranges strand
#          <Rle>           <IRanges>  <Rle>
#   [1]    chr12     8370000-8380000      *
#   [2]    chr13   92660000-92670000      *

dar_gr = gained_dar_gr
#dar_gr = lost_dar_gr
loop_anchor1 = loop_anchor1
loop_anchor2 = loop_anchor2


# ============================================================
# 1. Annotate loop anchors
# ============================================================

loop_anchor1_anno <- annotatePeak(
    loop_anchor1,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
)

loop_anchor2_anno <- annotatePeak(
    loop_anchor2,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
)

loop_anchor1_anno_gr <- loop_anchor1_anno@anno
loop_anchor2_anno_gr <- loop_anchor2_anno@anno


# ============================================================
# 2. Identify DARs overlapping each loop anchor
# ============================================================

anchor1_dar_overlap <- GrangeOverlaps(
    loop_anchor1,
    dar_gr,
    output = "ranges"
)

anchor2_dar_overlap <- GrangeOverlaps(
    loop_anchor2,
    dar_gr,
    output = "ranges"
)


# ============================================================
# 3. DAR at anchor 1 -> promoter at anchor 2
# ============================================================

idx_anchor1_dar <- (
    !is.na(anchor1_dar_overlap$grB) &
    is.na(anchor2_dar_overlap$grB) &
    grepl("Promoter", loop_anchor2_anno_gr$annotation)
)

anchor1_dar_gr <- anchor1_dar_overlap[
    idx_anchor1_dar,
]

anchor2_promoter_gr <- loop_anchor2_anno_gr[
    idx_anchor1_dar
]


# ============================================================
# 4. DAR at anchor 2 -> promoter at anchor 1
# ============================================================

idx_anchor2_dar <- (
    !is.na(anchor2_dar_overlap$grB) &
    is.na(anchor1_dar_overlap$grB) &
    grepl("Promoter", loop_anchor1_anno_gr$annotation)
)

anchor2_dar_gr <- anchor2_dar_overlap[
    idx_anchor2_dar,
]

anchor1_promoter_gr <- loop_anchor1_anno_gr[
    idx_anchor2_dar
]


# ============================================================
# 5. Collect loop-linked promoter target genes
# ============================================================

dar_loop_promoter_gr <- c(
    anchor2_promoter_gr,
    anchor1_promoter_gr
)

loop_target_genes <- unique(
    na.omit(dar_loop_promoter_gr$SYMBOL)
)

length(loop_target_genes)



# ============================================================
# 6. Generate DAR-loop-promoter mapping tables
# ============================================================
# Merge paired DAR-overlapping and promoter-overlapping anchors
# into a single table while retaining their one-to-one loop
# relationships.

merge_loop_anchors <- function(
    anchor1,
    anchor2,
    prefix1 = "anchor1_",
    prefix2 = "anchor2_"
) {

    # Convert GRanges objects to data frames
    df1 <- as.data.frame(anchor1)
    df2 <- as.data.frame(anchor2)

    # Check that paired anchors have the same number of rows
    if (nrow(df1) != nrow(df2)) {
        stop(
            paste0(
                "Anchor row numbers do not match: ",
                nrow(df1), " vs ", nrow(df2)
            )
        )
    }

    # Add prefixes to distinguish DAR and promoter anchor columns
    colnames(df1) <- paste0(prefix1, colnames(df1))
    colnames(df2) <- paste0(prefix2, colnames(df2))

    # Merge paired anchors
    combined_df <- cbind(df1, df2)

    # Clean column names generated by GrangeOverlaps()
    colnames(combined_df) <- gsub(
        "grA\\.|grB\\.",
        "",
        colnames(combined_df)
    )

    # Remove unnecessary columns
    combined_df[, c(
        "anchor1_width",
        "anchor1_strand",
        "anchor2_width",
        "anchor2_strand"
    )] <- NULL

    # Rename DAR identifier and sort by genomic position
    combined_df_sort <- combined_df %>%
        dplyr::rename(anchor1_DAR = anchor1_grB) %>%
        dplyr::arrange(anchor1_seqnames, anchor1_start)

    message(
        "Generated ",
        nrow(combined_df_sort),
        " DAR-loop-promoter links."
    )

    return(combined_df_sort)
}


# ============================================================
# 7. Merge both DAR-loop-promoter orientations
# ============================================================

# DAR at anchor 1 -> promoter at anchor 2
dar_anchor1_promoter_anchor2_df <- merge_loop_anchors(
    anchor1_dar_gr,
    anchor2_promoter_gr
)

# DAR at anchor 2 -> promoter at anchor 1
# Reorient the pair so that DAR is consistently represented
# as anchor 1 and promoter as anchor 2 in the output table.
dar_anchor2_promoter_anchor1_df <- merge_loop_anchors(
    anchor2_dar_gr,
    anchor1_promoter_gr
)

# Confirm identical table structure
stopifnot(
    identical(
        colnames(dar_anchor1_promoter_anchor2_df),
        colnames(dar_anchor2_promoter_anchor1_df)
    )
)

# Combine both loop orientations
dar_loop_target_df <- rbind(
    dar_anchor1_promoter_anchor2_df,
    dar_anchor2_promoter_anchor1_df
)

dim(dar_loop_target_df)
head(dar_loop_target_df)
dar_loop_target_df[1,]
#   anchor1_seqnames anchor1_start anchor1_end
# 1             chr1     114580000   114590000
#                                                                                                                                                                                                                                      anchor1_DAR
# 1 chr1\t114581177\t114582203\t1027\t+\t1027\tInterval_12719\t0.5961157\t0.0009429258\t0.03352546\tUp\tPromoter (<=1kb)\t1\t114567557\t114581615\t14059\t2\t10286\tENST00000369541.4\t0\tENSG00000116752\tBCAS2\tBCAS2 pre-mRNA processing factor
#   anchor2_seqnames anchor2_start anchor2_end anchor2_annotation anchor2_geneChr anchor2_geneStart anchor2_geneEnd anchor2_geneLength
# 1             chr1     114750000   114760000   Promoter (<=1kb)               1         114716913       114757983              41071
#   anchor2_geneStrand anchor2_geneId anchor2_transcriptId anchor2_distanceToTSS anchor2_ENSEMBL anchor2_SYMBOL                anchor2_GENENAME
# 1                  2           7812   ENST00000339438.10                     0 ENSG00000009307          CSDE1 cold shock domain containing E1




# ============================================================
# GRE-containing lost DAR target-gene analysis
# ============================================================
# This analysis identifies target genes associated with
# GRE-containing lost DARs through:
#
#   1. direct promoter overlap; and
#   2. chromatin loops connecting GRE-containing lost DARs
#      to promoter-overlapping anchors.
#
# Inputs:
#   GRE_regions_gr      : GRE-containing regions
#   lost_dar_gr         : annotated lost DARs
#   GRE_regions         : GRE-containing region table
#   dar_loop_target_df  : DAR-loop-promoter mapping table
# ============================================================


# ============================================================
# 1. Identify GRE-containing lost DARs
# ============================================================

GRE_dar_hits <- findOverlaps(
    GRE_regions_gr,
    lost_dar_gr,
    ignore.strand = TRUE
)

GRE_lost_dar_gr <- lost_dar_gr[
    unique(subjectHits(GRE_dar_hits))
]

length(GRE_lost_dar_gr)


# ============================================================
# 2. Identify promoter-associated target genes
# ============================================================

GRE_promoter_dar_gr <- GRE_lost_dar_gr[
    grepl("Promoter", GRE_lost_dar_gr$annotation)
]

GRE_promoter_target_genes <- unique(
    na.omit(GRE_promoter_dar_gr$SYMBOL)
)

length(GRE_promoter_target_genes)
GRE_promoter_target_genes


# ============================================================
# 3. Identify loop-linked target genes
# ============================================================
# GRE region starts were shifted by +1 during coordinate
# processing and therefore do not exactly match the original
# DAR coordinates. End coordinates are used to match GRE-
# containing DARs to the DAR-loop-promoter mapping table.

GRE_id <- as.character(end(GRE_regions_gr))

dar_id <- gsub("\t","_", dar_loop_target_df$anchor1_DAR)

GRE_loop_match <- Reduce(`|`,
    lapply(GRE_id,function(x) grepl(x,dar_id,fixed = TRUE)))

GRE_loop_target_df <- dar_loop_target_df[
    GRE_loop_match,
]

GRE_loop_target_genes <- unique(
    na.omit(GRE_loop_target_df$anchor2_SYMBOL)
)

length(GRE_loop_target_genes)
# NeuN+ lost DARs: 8



