library(dplyr)
library(readr)
library(survival)
library(survminer)
library(ggplot2)

dir.create("results/survival", recursive = TRUE, showWarnings = FALSE)

clinical <- read_tsv("data/processed/kirc_matched_clinical.tsv", show_col_types = FALSE)
clusters <- read_tsv("results/clustering/rppa_pca_clusters.tsv", show_col_types = FALSE)

surv_data <- clinical %>%
  select(patient_id, os_time, os_event, stage_clean, grade_clean) %>%
  inner_join(clusters, by = "patient_id") %>%
  filter(!is.na(os_time), !is.na(os_event), !is.na(cluster))

fit <- survfit(Surv(os_time, os_event) ~ cluster, data = surv_data)

png("results/survival/km_by_rppa_cluster.png", width = 900, height = 700)
print(
  ggsurvplot(
    fit,
    data = surv_data,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = FALSE,
    ggtheme = theme_minimal(),
    title = "Overall survival by RPPA-derived cluster in TCGA-KIRC"
  )
)
dev.off()

cluster_survival_summary <- surv_data %>%
  count(cluster, os_event)

write_tsv(cluster_survival_summary, "results/survival/cluster_survival_summary.tsv")

message("Survival analysis by cluster completed successfully.")