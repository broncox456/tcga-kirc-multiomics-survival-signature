library(dplyr)
library(readr)
library(tibble)

dir.create("results/integration", recursive = TRUE, showWarnings = FALSE)

clinical <- read_tsv("data/processed/kirc_clinical_clean.tsv", show_col_types = FALSE)
rna <- read_tsv("data/processed/kirc_rnaseq_patient_matrix.tsv", show_col_types = FALSE)
rppa <- read_tsv("data/processed/kirc_rppa.tsv", show_col_types = FALSE)

clinical_ids <- clinical %>%
  distinct(patient_id) %>%
  pull(patient_id)

rna_ids <- names(rna)[-1]

rppa_ids <- rppa %>%
  mutate(patient_id = substr(SampleID, 1, 12)) %>%
  distinct(patient_id) %>%
  pull(patient_id)

matched_ids <- Reduce(intersect, list(clinical_ids, rna_ids, rppa_ids))

matched_clinical <- clinical %>%
  filter(patient_id %in% matched_ids)

matched_rna <- rna %>%
  select(gene_id, all_of(matched_ids))

matched_rppa <- rppa %>%
  mutate(patient_id = substr(SampleID, 1, 12)) %>%
  filter(patient_id %in% matched_ids)

write_tsv(matched_clinical, "data/processed/kirc_matched_clinical.tsv")
write_tsv(matched_rna, "data/processed/kirc_matched_rnaseq.tsv")
write_tsv(matched_rppa, "data/processed/kirc_matched_rppa.tsv")

coverage <- tibble(
  layer = c("clinical", "rnaseq_patient_level", "rppa", "matched_all_three"),
  n = c(length(clinical_ids), length(rna_ids), length(unique(rppa_ids)), length(matched_ids))
)

write_tsv(coverage, "results/integration/multiomics_coverage.tsv")

message("Matched multi-omics cohort created successfully.")
message(paste("Matched patients across all three layers:", length(matched_ids)))