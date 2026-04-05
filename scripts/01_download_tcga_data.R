if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "dplyr",
  "readr",
  "stringr"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(readr)
library(stringr)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)

message("Downloading TCGA-KIRC clinical data...")
clinical <- GDCquery_clinic(project = "TCGA-KIRC", type = "clinical")

write_tsv(clinical, "data/metadata/tcga_kirc_clinical_raw.tsv")

message("Querying RNA-seq counts...")
query_rna <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

saveRDS(query_rna, "data/metadata/query_rna_tcga_kirc.rds")

message("Download script finished.")