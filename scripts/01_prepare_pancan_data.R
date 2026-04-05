library(dplyr)
library(readr)
library(readxl)
library(stringr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 1. CARGAR DATA
# -------------------------

rna <- read_tsv("data/raw/pancan_rnaseq.tsv")
rppa <- read_tsv("data/raw/pancan_rppa.txt")
clinical <- read_excel("data/raw/pancan_clinical.xlsx")

# -------------------------
# 2. FILTRAR KIRC
# -------------------------

clinical_kirc <- clinical %>%
  filter(type == "KIRC")

kirc_ids <- clinical_kirc$bcr_patient_barcode

# -------------------------
# 3. RNA (filtrar columnas KIRC)
# -------------------------

rna_kirc <- rna %>%
  select(Hugo_Symbol, one_of(kirc_ids))

# -------------------------
# 4. RPPA (filtrar KIRC)
# -------------------------

rppa_kirc <- rppa %>%
  select(Composite.Element.REF, one_of(kirc_ids))

# -------------------------
# 5. GUARDAR
# -------------------------

write_tsv(clinical_kirc, "data/processed/kirc_clinical.tsv")
write_tsv(rna_kirc, "data/processed/kirc_rnaseq.tsv")
write_tsv(rppa_kirc, "data/processed/kirc_rppa.tsv")

message("KIRC multi-omics datasets prepared successfully.")