library(dplyr)
library(readr)
library(stringr)
library(tibble)

dir.create("results/rnaseq", recursive = TRUE, showWarnings = FALSE)

rna <- read_tsv("data/processed/kirc_rnaseq.tsv", show_col_types = FALSE)

sample_cols <- names(rna)[-1]
patient_ids <- substr(sample_cols, 1, 12)

# escoger una sola muestra por paciente: la primera aparición
sample_map <- tibble(
  sample_col = sample_cols,
  patient_id = patient_ids
) %>%
  distinct(patient_id, .keep_all = TRUE)

selected_cols <- c("gene_id", sample_map$sample_col)

rna_patient <- rna %>%
  select(all_of(selected_cols))

# renombrar columnas a patient_id
new_names <- c("gene_id", sample_map$patient_id)
colnames(rna_patient) <- new_names

write_tsv(sample_map, "data/processed/kirc_rnaseq_sample_map.tsv")
write_tsv(rna_patient, "data/processed/kirc_rnaseq_patient_matrix.tsv")

rna_summary <- tibble(
  metric = c("n_genes", "n_unique_patients"),
  value = c(nrow(rna_patient), ncol(rna_patient) - 1)
)

write_tsv(rna_summary, "results/rnaseq/rnaseq_patient_summary.tsv")

message("RNA patient-level matrix created successfully.")