library(dplyr)
library(readr)
library(ggplot2)
library(tibble)

dir.create("results/clustering", recursive = TRUE, showWarnings = FALSE)

rppa <- read_tsv("data/processed/kirc_matched_rppa.tsv", show_col_types = FALSE)

rppa_matrix <- rppa %>%
  mutate(patient_id = substr(SampleID, 1, 12)) %>%
  select(-SampleID, -TumorType) %>%
  group_by(patient_id) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()

patient_ids <- rppa_matrix$patient_id

mat <- rppa_matrix %>%
  select(-patient_id) %>%
  as.data.frame()

rownames(mat) <- patient_ids

# quitar columnas totalmente NA o con varianza 0
mat <- mat[, colSums(is.na(mat)) < nrow(mat), drop = FALSE]
mat <- mat[, apply(mat, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]

# imputación simple por mediana
for (j in seq_len(ncol(mat))) {
  x <- mat[, j]
  x[is.na(x)] <- median(x, na.rm = TRUE)
  mat[, j] <- x
}

# escalar
mat_scaled <- scale(mat)

# PCA
pca <- prcomp(mat_scaled, center = TRUE, scale. = FALSE)

pca_df <- data.frame(
  patient_id = rownames(mat_scaled),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

# clustering kmeans simple
set.seed(123)
km <- kmeans(pca_scores[, c("PC1", "PC2")], centers = 2, nstart = 50)

pca_df$cluster <- as.factor(km$cluster)

write_tsv(pca_df, "results/clustering/rppa_pca_clusters.tsv")

explained_var <- summary(pca)$importance[2, 1:2] * 100

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "TCGA-KIRC RPPA PCA with k-means clusters",
    x = paste0("PC1 (", round(explained_var[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(
  filename = "results/clustering/rppa_pca_clusters.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

cluster_sizes <- pca_df %>%
  count(cluster, sort = TRUE)

write_tsv(cluster_sizes, "results/clustering/rppa_cluster_sizes.tsv")

message("RPPA PCA and clustering completed successfully.")
message(paste("Patients clustered:", nrow(pca_df)))