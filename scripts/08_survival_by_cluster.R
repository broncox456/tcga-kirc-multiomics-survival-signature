library(readr)
library(dplyr)
library(survival)

dir.create("results/survival", recursive = TRUE, showWarnings = FALSE)

clusters <- read_tsv("results/clustering/rppa_pca_clusters.tsv", show_col_types = FALSE)
clinical <- read_tsv("data/processed/kirc_clinical_clean.tsv", show_col_types = FALSE)

# unir clusters con supervivencia clínica
surv_data <- clusters %>%
  inner_join(
    clinical %>% select(patient_id, os_time, os_event, stage_clean, grade_clean),
    by = "patient_id"
  ) %>%
  filter(!is.na(os_time), !is.na(os_event), !is.na(cluster)) %>%
  mutate(
    os_time = as.numeric(os_time),
    os_event = as.numeric(os_event),
    cluster = as.factor(cluster)
  )

fit <- survfit(Surv(os_time, os_event) ~ cluster, data = surv_data)

png("results/survival/kaplan_meier_by_cluster.png", width = 1000, height = 800)

plot(
  fit,
  col = c("red", "blue", "darkgreen"),
  lwd = 2,
  xlab = "Time (days)",
  ylab = "Survival probability",
  main = "Kaplan-Meier Survival by RPPA Clusters"
)

legend(
  "bottomleft",
  legend = levels(surv_data$cluster),
  col = c("red", "blue", "darkgreen"),
  lwd = 2,
  bty = "n"
)

dev.off()

surv_summary <- surv_data %>%
  count(cluster, os_event)

write_tsv(surv_summary, "results/survival/survival_cluster_counts.tsv")

message("Kaplan-Meier analysis completed successfully.")
message(paste("Patients used in survival analysis:", nrow(surv_data)))