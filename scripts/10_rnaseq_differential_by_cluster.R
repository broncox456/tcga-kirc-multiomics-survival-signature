library(readr)
library(dplyr)
library(tibble)

dir.create("results/differential", recursive = TRUE, showWarnings = FALSE)

clusters <- read_tsv("results/clustering/rppa_pca_clusters.tsv", show_col_types = FALSE) %>%
  select(patient_id, cluster) %>%
  mutate(cluster = as.factor(cluster))

rna <- read_tsv("data/processed/kirc_matched_rnaseq.tsv", show_col_types = FALSE)

# matriz genes x pacientes
gene_ids <- rna$gene_id
rna_mat <- as.data.frame(rna[, -1])
rownames(rna_mat) <- gene_ids

# pacientes comunes
common_ids <- intersect(colnames(rna_mat), clusters$patient_id)
rna_mat <- rna_mat[, common_ids, drop = FALSE]
clusters <- clusters %>% filter(patient_id %in% common_ids)

# ordenar igual
clusters <- clusters[match(common_ids, clusters$patient_id), ]
stopifnot(all(clusters$patient_id == colnames(rna_mat)))

# log2
rna_mat[] <- lapply(rna_mat, as.numeric)
rna_log <- log2(as.matrix(rna_mat) + 1)

# limpiar nombres de genes
gene_symbols <- rownames(rna_log)
gene_symbols <- sub("\\|.*$", "", gene_symbols)
gene_symbols <- sub("\\..*$", "", gene_symbols)

group1_ids <- clusters$patient_id[clusters$cluster == "1"]
group2_ids <- clusters$patient_id[clusters$cluster == "2"]

results <- lapply(seq_len(nrow(rna_log)), function(i) {
  x1 <- as.numeric(rna_log[i, group1_ids])
  x2 <- as.numeric(rna_log[i, group2_ids])

  test <- tryCatch(
    wilcox.test(x1, x2),
    error = function(e) NULL
  )

  tibble(
    gene_id = rownames(rna_log)[i],
    gene_symbol = gene_symbols[i],
    median_cluster1 = median(x1, na.rm = TRUE),
    median_cluster2 = median(x2, na.rm = TRUE),
    diff = median_cluster1 - median_cluster2,
    p_value = if (!is.null(test)) test$p.value else NA_real_
  )
})

rna_diff <- bind_rows(results) %>%
  filter(!is.na(p_value), !is.na(gene_symbol), gene_symbol != "") %>%
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_value)

write_tsv(rna_diff, "results/differential/rnaseq_differential.tsv")

top_up_cluster1 <- rna_diff %>%
  arrange(desc(diff), p_value) %>%
  slice_head(n = 30)

top_up_cluster2 <- rna_diff %>%
  arrange(diff, p_value) %>%
  slice_head(n = 30)

write_tsv(top_up_cluster1, "results/differential/top_cluster1_genes.tsv")
write_tsv(top_up_cluster2, "results/differential/top_cluster2_genes.tsv")

message("RNA differential analysis completed successfully.")
message(paste("Genes tested:", nrow(rna_diff)))