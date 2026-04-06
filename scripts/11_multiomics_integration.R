library(readr)
library(dplyr)
library(tibble)

dir.create("results/integration", recursive = TRUE, showWarnings = FALSE)

rppa_diff <- read_tsv("results/differential/rppa_differential.tsv", show_col_types = FALSE)
rna_diff  <- read_tsv("results/differential/rnaseq_differential.tsv", show_col_types = FALSE)

# mapeo manual curado de proteínas RPPA -> gen canónico
protein_gene_map <- tribble(
  ~protein,         ~gene_symbol, ~pathway_theme,
  "PAI1",           "SERPINE1",   "ECM / invasion",
  "BAK",            "BAK1",       "Apoptosis",
  "RAD51",          "RAD51",      "DNA repair",
  "ASNS",           "ASNS",       "Metabolism",
  "IGFBP2",         "IGFBP2",     "Growth signaling",
  "P53",            "TP53",       "Cell cycle / stress",
  "CYCLINB1",       "CCNB1",      "Cell cycle",
  "STATHMIN",       "STMN1",      "Proliferation / mitosis",
  "CYCLINE1",       "CCNE1",      "Cell cycle",
  "BCLXL",          "BCL2L1",     "Apoptosis",
  "TFRC",           "TFRC",       "Iron metabolism / growth",
  "SMAC",           "DIABLO",     "Apoptosis",
  "GATA3",          "GATA3",      "Transcriptional regulation",
  "P70S6K_pT389",   "RPS6KB1",    "mTOR signaling",
  "CHK2_pT68",      "CHEK2",      "DNA damage response",
  "BRCA2",          "BRCA2",      "DNA repair",
  "P21",            "CDKN1A",     "Cell cycle arrest"
)

integration <- protein_gene_map %>%
  inner_join(rppa_diff %>% select(protein, median_cluster1, median_cluster2, diff, p_value, fdr),
             by = "protein") %>%
  rename(
    protein_median_cluster1 = median_cluster1,
    protein_median_cluster2 = median_cluster2,
    protein_diff = diff,
    protein_p_value = p_value,
    protein_fdr = fdr
  ) %>%
  inner_join(rna_diff %>% select(gene_symbol, median_cluster1, median_cluster2, diff, p_value, fdr),
             by = "gene_symbol") %>%
  rename(
    rna_median_cluster1 = median_cluster1,
    rna_median_cluster2 = median_cluster2,
    rna_diff = diff,
    rna_p_value = p_value,
    rna_fdr = fdr
  ) %>%
  mutate(
    concordant_direction = case_when(
      protein_diff > 0 & rna_diff > 0 ~ "Up in Cluster 1 at both layers",
      protein_diff < 0 & rna_diff < 0 ~ "Up in Cluster 2 at both layers",
      TRUE ~ "Discordant"
    )
  ) %>%
  arrange(concordant_direction, protein_p_value, rna_p_value)

write_tsv(integration, "results/integration/multiomics_concordance.tsv")

top_concordant_cluster1 <- integration %>%
  filter(concordant_direction == "Up in Cluster 1 at both layers") %>%
  arrange(protein_p_value, rna_p_value)

write_tsv(top_concordant_cluster1, "results/integration/top_concordant_cluster1.tsv")

message("Multi-omics integration completed successfully.")
message(paste("Integrated markers:", nrow(integration)))