library(dplyr)
library(readr)
library(readxl)
library(stringr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 1. CARGAR DATA
# -------------------------

rna <- read_tsv("data/raw/pancan_rnaseq.tsv", show_col_types = FALSE)
rppa <- read_tsv("data/raw/pancan_rppa.txt", show_col_types = FALSE)
clinical <- read_excel("data/raw/pancan_clinical.xlsx")

# -------------------------
# 2. FILTRAR KIRC EN CLINICA
# -------------------------

clinical_kirc <- clinical %>%
  filter(type == "KIRC")

kirc_ids <- unique(clinical_kirc$bcr_patient_barcode)
kirc_ids <- kirc_ids[!is.na(kirc_ids)]

# -------------------------
# 3. PREPARAR RNA
# El archivo RNA tiene gene_id y columnas con barcodes largos
# -------------------------

rna_sample_cols <- names(rna)[-1]
rna_patient_ids <- substr(rna_sample_cols, 1, 12)

rna_keep_cols <- c(
  "gene_id",
  rna_sample_cols[rna_patient_ids %in% kirc_ids]
)

rna_kirc <- rna %>%
  select(all_of(rna_keep_cols))

# -------------------------
# 4. PREPARAR RPPA
# El archivo RPPA viene en formato filas = muestras
# -------------------------

rppa_kirc <- rppa %>%
  filter(TumorType == "KIRC") %>%
  mutate(
    patient_id = substr(SampleID, 1, 12)
  ) %>%
  filter(patient_id %in% kirc_ids)

# -------------------------
# 5. GUARDAR
# -------------------------

write_tsv(clinical_kirc, "data/processed/kirc_clinical.tsv")
write_tsv(rna_kirc, "data/processed/kirc_rnaseq.tsv")
write_tsv(rppa_kirc, "data/processed/kirc_rppa.tsv")

message("KIRC multi-omics datasets prepared successfully.")
message(paste("Clinical KIRC patients:", nrow(clinical_kirc)))
message(paste("RNA KIRC columns:", ncol(rna_kirc) - 1))
message(paste("RPPA KIRC samples:", nrow(rppa_kirc)))