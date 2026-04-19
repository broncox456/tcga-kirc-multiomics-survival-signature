library(SummarizedExperiment)
library(dplyr)
library(readr)
library(stringr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/rnaseq", recursive = TRUE, showWarnings = FALSE)

rna_se <- readRDS("data/raw/tcga_kirc_rnaseq_se.rds")

# matriz de counts
count_matrix <- assay(rna_se)

# metadatos de muestras
sample_meta <- as.data.frame(colData(rna_se))

sample_meta <- sample_meta %>%
  mutate(
    sample_barcode = substr(colnames(count_matrix), 1, 16),
    patient_id = substr(colnames(count_matrix), 1, 12)
  )

# guardar sample metadata
write_tsv(sample_meta, "data/processed/tcga_kirc_rnaseq_sample_metadata.tsv")

# filtrar genes de baja expresión simple
keep_genes <- rowSums(count_matrix >= 10) >= 10
count_matrix_filtered <- count_matrix[keep_genes, ]

saveRDS(count_matrix_filtered, "data/processed/tcga_kirc_counts_filtered.rds")

# tabla resumen
rna_summary <- tibble::tibble(
  metric = c("n_genes_raw", "n_genes_filtered", "n_samples"),
  value = c(nrow(count_matrix), nrow(count_matrix_filtered), ncol(count_matrix_filtered))
)

write_tsv(rna_summary, "results/rnaseq/rnaseq_summary.tsv")

message("RNA-seq matrix prepared successfully.")