library(readr)
library(dplyr)
library(ggplot2)
library(tibble)

dir.create("results/clustering", recursive = TRUE, showWarnings = FALSE)

# Cargar RPPA matched
rppa <- read_tsv("data/processed/kirc_matched_rppa.tsv", show_col_types = FALSE)

# Construir matriz paciente × proteínas
rppa_patient <- rppa %>%
  mutate(patient_id = substr(SampleID, 1, 12)) %>%
  select(-SampleID, -TumorType) %>%
  group_by(patient_id) %>%
  summarise(across(everything(), ~ mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

rppa_mat <- as.data.frame(rppa_patient[, -1])
rownames(rppa_mat) <- rppa_patient$patient_id

# Quitar columnas totalmente NA
keep_cols <- colSums(is.na(rppa_mat)) < nrow(rppa_mat)
rppa_mat <- rppa_mat[, keep_cols, drop = FALSE]

# Imputación simple por mediana
for (j in seq_len(ncol(rppa_mat))) {
  na_idx <- is.na(rppa_mat[, j])
  if (any(na_idx)) {
    med_val <- median(rppa_mat[, j], na.rm = TRUE)
    rppa_mat[na_idx, j] <- med_val
  }
}

# Quitar proteínas sin varianza
rppa_var <- apply(rppa_mat, 2, var, na.rm = TRUE)
rppa_mat <- rppa_mat[, rppa_var > 0, drop = FALSE]

# PCA
pca_res <- prcomp(rppa_mat, scale. = TRUE)

pca_scores <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("patient_id")

# Clustering con k = 2
set.seed(123)
km <- kmeans(pca_scores[, c("PC1", "PC2")], centers = 2, nstart = 50)

pca_scores$cluster <- as.factor(km$cluster)

# Guardar tabla de clusters
write_tsv(pca_scores, "results/clustering/rppa_pca_clusters.tsv")

cluster_sizes <- pca_scores %>%
  count(cluster)

write_tsv(cluster_sizes, "results/clustering/rppa_cluster_sizes.tsv")

# Gráfico PCA
p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    title = "RPPA PCA Clustering in TCGA-KIRC",
    x = "PC1",
    y = "PC2",
    color = "Cluster"
  ) +
  theme_minimal()

ggsave("results/clustering/rppa_pca_clusters.png", p, width = 8, height = 6, dpi = 300)

message("RPPA PCA and clustering completed successfully.")
message(paste("Patients clustered:", nrow(pca_scores)))
print(cluster_sizes)