library(readr)
library(dplyr)
library(tibble)

dir.create("results/differential", recursive = TRUE, showWarnings = FALSE)

clusters <- read_tsv("results/clustering/rppa_pca_clusters.tsv", show_col_types = FALSE)
rppa <- read_tsv("data/processed/kirc_matched_rppa.tsv", show_col_types = FALSE)

# preparar matriz paciente
rppa_patient <- rppa %>%
  mutate(patient_id = substr(SampleID, 1, 12)) %>%
  select(-SampleID, -TumorType) %>%
  group_by(patient_id) %>%
  summarise(across(everything(), ~ mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

# unir con clusters
data <- rppa_patient %>%
  inner_join(clusters, by = "patient_id")

# asegurar cluster numérico/entero simple
data <- data %>%
  mutate(cluster = as.integer(as.character(cluster)))

# proteínas reales solamente
proteins <- setdiff(names(data), c("patient_id", "cluster", "PC1", "PC2"))

# separar grupos
group1 <- data %>% filter(cluster == 1)
group2 <- data %>% filter(cluster == 2)

results <- lapply(proteins, function(p) {
  x1 <- suppressWarnings(as.numeric(group1[[p]]))
  x2 <- suppressWarnings(as.numeric(group2[[p]]))

  test <- tryCatch(
    wilcox.test(x1, x2),
    error = function(e) NULL
  )

  tibble(
    protein = p,
    median_cluster1 = median(x1, na.rm = TRUE),
    median_cluster2 = median(x2, na.rm = TRUE),
    diff = median_cluster1 - median_cluster2,
    p_value = if (!is.null(test)) test$p.value else NA_real_
  )
})

results_df <- bind_rows(results) %>%
  filter(!is.na(p_value)) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write_tsv(results_df, "results/differential/rppa_differential.tsv")

top_up_cluster1 <- results_df %>%
  filter(diff > 0) %>%
  arrange(p_value, desc(diff)) %>%
  slice_head(n = 20)

top_up_cluster2 <- results_df %>%
  filter(diff < 0) %>%
  arrange(p_value, diff) %>%
  slice_head(n = 20)

write_tsv(top_up_cluster1, "results/differential/top_cluster1_proteins.tsv")
write_tsv(top_up_cluster2, "results/differential/top_cluster2_proteins.tsv")

message("Differential protein analysis completed successfully.")
message(paste("Proteins tested:", nrow(results_df)))