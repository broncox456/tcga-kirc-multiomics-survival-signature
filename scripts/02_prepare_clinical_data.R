library(dplyr)
library(readr)
library(stringr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/clinical", recursive = TRUE, showWarnings = FALSE)

clinical <- read_tsv("data/metadata/tcga_kirc_clinical_raw.tsv", show_col_types = FALSE)

# inspección básica
write_tsv(
  tibble::tibble(column_name = names(clinical)),
  "results/clinical/clinical_columns_inventory.tsv"
)

# seleccionar variables clínicamente útiles
clinical_clean <- clinical %>%
  mutate(
    patient_id = submitter_id,
    age_at_diagnosis = as.numeric(age_at_diagnosis),
    days_to_death = suppressWarnings(as.numeric(days_to_death)),
    days_to_last_follow_up = suppressWarnings(as.numeric(days_to_last_follow_up)),
    os_time = dplyr::coalesce(days_to_death, days_to_last_follow_up),
    os_event = case_when(
      vital_status == "Dead" ~ 1,
      vital_status == "Alive" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  select(
    patient_id,
    age_at_diagnosis,
    gender,
    race,
    ajcc_pathologic_stage,
    ajcc_pathologic_t,
    ajcc_pathologic_n,
    ajcc_pathologic_m,
    tumor_grade,
    vital_status,
    days_to_death,
    days_to_last_follow_up,
    os_time,
    os_event
  )

write_tsv(clinical_clean, "data/processed/tcga_kirc_clinical_clean.tsv")

clinical_summary <- clinical_clean %>%
  summarise(
    n_patients = n(),
    n_with_os_time = sum(!is.na(os_time)),
    n_events = sum(os_event == 1, na.rm = TRUE),
    median_age = median(age_at_diagnosis, na.rm = TRUE)
  )

write_tsv(clinical_summary, "results/clinical/clinical_summary.tsv")

stage_distribution <- clinical_clean %>%
  count(ajcc_pathologic_stage, sort = TRUE)

write_tsv(stage_distribution, "results/clinical/stage_distribution.tsv")

message("Clinical data cleaned successfully.")