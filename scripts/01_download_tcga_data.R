if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_packages <- c("dplyr", "readr", "stringr", "tibble")
bioc_packages <- c("TCGAbiolinks", "SummarizedExperiment")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)

message("Step 1: Downloading TCGA-KIRC clinical data...")
clinical <- GDCquery_clinic(project = "TCGA-KIRC", type = "clinical")
write_tsv(clinical, "data/metadata/tcga_kirc_clinical_raw.tsv")

message("Step 2: Querying RNA-seq counts...")
query_rna <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

saveRDS(query_rna, "data/metadata/query_rna_tcga_kirc.rds")

message("Step 3: Downloading RNA-seq counts from GDC...")
GDCdownload(query_rna)

message("Step 4: Preparing SummarizedExperiment object...")
rna_se <- GDCprepare(query_rna)

saveRDS(rna_se, "data/raw/tcga_kirc_rnaseq_se.rds")

message("Done: clinical and RNA-seq data saved successfully.")