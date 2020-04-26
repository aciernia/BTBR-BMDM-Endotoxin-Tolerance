Code for differential analysis of ATAC-Seq peaks for strain and LPS treatment effects.

ATAC_DiffBindanalysis_4_27_19.R code for identify consensus peaks and differential peaks using DiffBind: https://bioconductor.org/packages/release/bioc/html/DiffBind.html

Uses library depth and RLE normalization, and EdgeR/limma-voom to identify differential peaks based on normalized ATAC-signal using a statistical model for strain and LPS treatment with the correlation for the within subjects variation removed. Also has additional QC measures.

OverlapComparisons_1_5_2020.R DAR parsing and permutation testing for overlapping region lists between strain conditions.

GenomicElementEnrichment.R Annotation of DARs using ChIPSeeker and testing for enrichment of DARs within genomic elements compared to all consensus background peaks using a Fisher's Exact test with FDR correction.