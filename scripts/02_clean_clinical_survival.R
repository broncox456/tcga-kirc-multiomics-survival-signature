library(dplyr)
library(readr)
library(tibble)

dir.create("results/clinical", recursive = TRUE, showWarnings = FALSE)

clinical <- read_tsv("data/processed/kirc_clinical.tsv", show_col_types = FALSE)

clinical_clean <- clinical %>%
  mutate(
    patient_id = bcr_patient_barcode,
    age_clean = suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis)),
    os_time = suppressWarnings(as.numeric(OS.time)),
    os_event = suppressWarnings(as.numeric(OS)),
    death_days_to_clean = suppressWarnings(as.numeric(death_days_to)),
    last_contact_days_to_clean = suppressWarnings(as.numeric(last_contact_days_to)),
    stage_clean = as.character(clinical_stage),
    grade_clean = as.character(histological_grade),
    vital_status_clean = as.character(vital_status)
  ) %>%
  distinct(patient_id, .keep_all = TRUE)

write_tsv(clinical_clean, "data/processed/kirc_clinical_clean.tsv")

clinical_summary <- clinical_clean %>%
  summarise(
    n_patients = n(),
    n_with_os_time = sum(!is.na(os_time)),
    n_events = sum(os_event == 1, na.rm = TRUE),
    median_age = median(age_clean, na.rm = TRUE)
  )

write_tsv(clinical_summary, "results/clinical/clinical_summary.tsv")

stage_distribution <- clinical_clean %>%
  count(stage_clean, sort = TRUE)

write_tsv(stage_distribution, "results/clinical/stage_distribution.tsv")

grade_distribution <- clinical_clean %>%
  count(grade_clean, sort = TRUE)

write_tsv(grade_distribution, "results/clinical/grade_distribution.tsv")

message("Clinical survival dataset cleaned successfully.")
message(paste("Patients:", nrow(clinical_clean)))
message(paste("Patients with OS.time:", sum(!is.na(clinical_clean$os_time))))
message(paste("Events:", sum(clinical_clean$os_event == 1, na.rm = TRUE)))