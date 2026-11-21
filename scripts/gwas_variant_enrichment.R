# ============================================================
# PTSD risk-variant enrichment in DARs
# ============================================================
# This analysis tests whether PTSD-associated risk variants
# overlap DARs more frequently than expected by chance.
#
# DARs are randomized within the same chromosome while preserving
# their genomic widths. ENCODE blacklist regions are excluded.
# The observed number of DAR-risk variant overlaps is compared
# with the distribution obtained from 1,000 permutations.
#
# Inputs:
#   dar_gr       : DARs used for enrichment analysis
#   risk_snp_gr  : PTSD risk variants in hg38 coordinates
# ============================================================

library(GenomicRanges)
library(GenomeInfoDb)
library(regioneR)
library(pbapply)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)



dar_gr
# GRanges object with 6415 ranges and 2 metadata columns:
#        seqnames              ranges strand |    Length          Geneid
#           <Rle>           <IRanges>  <Rle> | <integer>     <character>
#   3845     chr8   12192993-12195168      + |      2176 Interval_255163
#   1640    chr21     8392653-8394462      + |      1810 Interval_159625
#   3529     chr7   56373408-56373893      + |       486 Interval_243471
#   2621     chr4 189233835-189234182      + |       348 Interval_203533
risk_snp_gr
# GRanges object with 12901 ranges and 2 metadata columns:
#           seqnames    ranges strand |         SNP         P
#              <Rle> <IRanges>  <Rle> | <character> <numeric>
#       [1]     chr2  22207923      * |   rs1554266 4.357e-18
#       [2]     chr2  22208074      * |   rs1554267 7.963e-18


# ============================================================
# 1. Define genomic background
# ============================================================

chrom_sizes <- seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)

keep_chr <- paste0("chr",c(1:22, "X", "Y"))

chrom_sizes <- chrom_sizes[keep_chr]

genome_gr <- GRanges(
    seqnames = names(chrom_sizes),
    ranges = IRanges(
        start = 1,
        end = as.numeric(chrom_sizes)
    )
)


# ============================================================
# 2. Load hg38 blacklist regions
# ============================================================

blacklist_file <- "./ENCFF356LFX_hg38_blacklist_v3.bed"

blacklist_df <- read.table(
    blacklist_file,
    header = FALSE,
    stringsAsFactors = FALSE
)

colnames(blacklist_df) <- c(
    "seqnames",
    "start",
    "end"
)

blacklist_df$seqnames <- paste0(
    "chr",
    blacklist_df$seqnames
)


# ============================================================
# 3. Calculate observed DAR-risk variant overlaps
# ============================================================

ov_index <-  findOverlaps(
        dar_gr,
        risk_snp_gr,
        ignore.strand = TRUE)

#observed_overlaps <- length(unique(subjectHits(ov_index)))
observed_overlaps <- length(ov_index)
observed_overlaps
#[1] 62 NP
#[1] 23 NN


# ============================================================
# 4. Permute DAR locations
# ============================================================

set.seed(123)

n_permutations <- 1000

random_overlaps <- pbreplicate(
    n_permutations,
    {
        random_dar_gr <- randomizeRegions(
            dar_gr,
            genome = genome_gr,
            mask = blacklist_df,
            per.chromosome = TRUE,
            allow.overlaps = FALSE
        )

        length(
            findOverlaps(
                random_dar_gr,
                risk_snp_gr,
                ignore.strand = TRUE
            )
        )
    }
)


# ============================================================
# 5. Calculate empirical enrichment P value
# ============================================================

p_empirical <- mean(random_overlaps >= observed_overlaps)
p_empirical






# ============================================================
# Link distal PTSD risk variants to promoter target genes
# ============================================================
# Identifies genes whose promoter-overlapping loop anchors are
# connected to distal anchors containing PTSD risk variants.
#
# Inputs:
#   loop_anchor1 : first loop anchors
#   loop_anchor2 : second loop anchors
#   risk_snp_gr  : PTSD risk variants
#
# Outputs:
#   mapping      : SNP-loop-promoter mapping table
#   gene_summary : target genes and linked SNP counts
# ============================================================

link_distal_SNP_to_promoter <- function(loop_anchor1, loop_anchor2, sig_snp) {

    # Annotate loop anchors
    loop_anchor1_anno <- annotatePeak(loop_anchor1, tssRegion = c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db", verbose = FALSE)
    loop_anchor2_anno <- annotatePeak(loop_anchor2, tssRegion = c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db", verbose = FALSE)

    # Identify PTSD risk variants overlapping loop anchors
    anchor1_snp_overlap <- GrangeOverlaps(loop_anchor1_anno@anno, sig_snp, output = "ranges")
    anchor2_snp_overlap <- GrangeOverlaps(loop_anchor2_anno@anno, sig_snp, output = "ranges")

    # Convert to data frames
    anchor1_df <- data.frame(data.frame(granges(anchor1_snp_overlap$grA.X))[, 1:3], anchor1_snp_overlap[, -1])
    anchor2_df <- data.frame(data.frame(granges(anchor2_snp_overlap$grA.X))[, 1:3], anchor2_snp_overlap[, -1])

    colnames(anchor1_df) <- paste0(colnames(anchor1_df), "_1")
    colnames(anchor2_df) <- paste0(colnames(anchor2_df), "_2")

    loop_snp_df <- cbind(anchor1_df, anchor2_df)

    # Retain loops containing PTSD risk variants
    loop_snp_df <- loop_snp_df %>%
        dplyr::filter(!(is.na(grB_1) & is.na(grB_2)))

    # Identify distal SNP-to-promoter interactions
    snp_loop_promoter_df <- loop_snp_df %>%
        dplyr::filter(
            (grepl("Promoter", grA.annotation_1) & !is.na(grB_2)) |
            (grepl("Promoter", grA.annotation_2) & !is.na(grB_1))
        ) %>%
        dplyr::select(seqnames_1, start_1, end_1, grA.SYMBOL_1, grA.annotation_1, grB_1,
                      seqnames_2, start_2, end_2, grA.SYMBOL_2, grA.annotation_2, grB_2)

    # Assign promoter target genes and count linked SNPs
    distal_snp_promoter_df <- snp_loop_promoter_df %>%
        dplyr::mutate(
            gene_keep1 = ifelse(grepl("Promoter", grA.annotation_1) & !is.na(grB_2), grA.SYMBOL_1, NA),
            snp_count1 = ifelse(grepl("Promoter", grA.annotation_1) & !is.na(grB_2), stringr::str_count(grB_2, "rs"), 0),
            gene_keep2 = ifelse(grepl("Promoter", grA.annotation_2) & !is.na(grB_1), grA.SYMBOL_2, NA),
            snp_count2 = ifelse(grepl("Promoter", grA.annotation_2) & !is.na(grB_1), stringr::str_count(grB_1, "rs"), 0)
        )

    # Summarize linked risk variants by target gene
    target_gene_summary <- dplyr::bind_rows(
        distal_snp_promoter_df %>% dplyr::select(gene = gene_keep1, snp_count = snp_count1),
        distal_snp_promoter_df %>% dplyr::select(gene = gene_keep2, snp_count = snp_count2)
    ) %>%
        dplyr::filter(!is.na(gene)) %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(total_snp_count = sum(snp_count), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(total_snp_count))

    return(list(mapping = distal_snp_promoter_df, gene_summary = target_gene_summary))
}


distal_SNP_to_promoter_list = link_distal_SNP_to_promoter(loop_anchor1, loop_anchor2, risk_snp_gr)




# ============================================================
# PTSD risk-variant enrichment in loop-mediated interactions
# ============================================================
# Tests whether observed chromatin loops connect promoter regions
# to more distal PTSD risk variants than expected by chance.
#
# Loop spans are randomized within chromosomes while preserving
# their genomic lengths. Randomized loop anchors are reconstructed
# at 10-kb resolution, and distal SNP-to-promoter links are
# quantified using link_distalSNP_to_promoter().
#
# Inputs:
#   loop_anchor1 : first loop anchors
#   loop_anchor2 : second loop anchors
#   risk_snp_gr  : PTSD risk variants
# ============================================================


# ============================================================
# 1. Construct loop spans
# ============================================================

starts <- pmin(end(loop_anchor1), start(loop_anchor2))
ends <- pmax(end(loop_anchor1), start(loop_anchor2))

loop_gr <- GRanges(
    seqnames = seqnames(loop_anchor1),
    ranges = IRanges(start = starts, end = ends)
)

length(loop_gr)


# ============================================================
# 2. Calculate observed SNP-loop-promoter links
# ============================================================

#observed_result <- link_distal_SNP_to_promoter(loop_anchor1, loop_anchor2, risk_snp_gr)

n_overlap_real <- sum(
    distal_SNP_to_promoter_list$gene_summary$total_snp_count,
    na.rm = TRUE
)

n_overlap_real

# NeuN+: 896
# NeuN-: 495


# ============================================================
# 3. Randomize loop locations
# ============================================================

loop_res <- 10000
n_permutations <- 1000

set.seed(123)

n_overlap_random <- pbreplicate(n_permutations, {

    random_loop_gr <- randomizeRegions(
        loop_gr,
        genome = genome_gr,
        mask = blacklist_df,
        per.chromosome = TRUE,
        allow.overlaps = TRUE
    )

    random_loop_anchor1 <- GRanges(
        seqnames = seqnames(random_loop_gr),
        ranges = IRanges(
            start = start(random_loop_gr) - loop_res,
            end = start(random_loop_gr)
        )
    )

    random_loop_anchor2 <- GRanges(
        seqnames = seqnames(random_loop_gr),
        ranges = IRanges(
            start = end(random_loop_gr),
            end = end(random_loop_gr) + loop_res
        )
    )

    random_result <- link_distal_SNP_to_promoter(
        random_loop_anchor1,
        random_loop_anchor2,
        risk_snp_gr
    )

    sum(random_result$gene_summary$total_snp_count, na.rm = TRUE)
})


# ============================================================
# 4. Calculate empirical P value
# ============================================================

p_val <- mean(n_overlap_random >= n_overlap_real)

p_val


