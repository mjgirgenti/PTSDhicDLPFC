# Code for “Multiscale chromatin alterations linked to transcriptional dysregulation in human PTSD”

This repository contains the main analysis code associated with the study.

The study investigates cell-type-specific alterations in chromatin accessibility and three-dimensional (3D) genome organization in the dorsolateral prefrontal cortex (DLPFC) of individuals with posttraumatic stress disorder (PTSD) and unaffected controls.

Bulk ATAC-seq and Hi-C data from neuronal (NeuN+) and non-neuronal (NeuN−) nuclei were integrated with transcriptomic and genetic data to characterize differential chromatin accessibility, chromatin looping, topologically associating domain (TAD) organization, and their relationships with PTSD-associated genetic variation and gene expression.

## Code

The repository contains the following R scripts:

- `dar_analysis.R`  
  Differential chromatin accessibility analysis of ATAC-seq data.

- `dar_loop_integration.R`  
  Integration of differentially accessible regions (DARs) with differential chromatin loops.

- `dar_target_gene_gre_analysis.R`  
  Analysis of DAR-associated target genes and glucocorticoid response element (GRE)-related regulatory features.

- `gwas_variant_enrichment.R`  
  Enrichment analysis of PTSD GWAS risk variants in chromatin accessibility and 3D genome features.

- `differential_loop.R`  
  Differential chromatin loop analysis of Hi-C data.

- `strengthened_loop_tad_context.R`  
  Analysis of PTSD-associated strengthened chromatin loops in the context of TAD organization.

- `differential_tad_boundary.R`  
  Differential TAD boundary analysis of Hi-C data.

- `tad_boundary_feature_enrichment.R`  
  Enrichment analysis of regulatory and molecular features near differential TAD boundaries.

- `tad_boundary_hotspot_analysis.R`  
  Analysis of genomic regions containing clusters of differential TAD boundaries.

- `tad_dar_correlation.R`  
  Analysis of relationships between TAD boundary alterations and differential chromatin accessibility.


