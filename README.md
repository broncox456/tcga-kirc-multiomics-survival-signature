# 🧬 TCGA-KIRC Multi-Omics Survival Stratification

![Kaplan-Meier Survival](results/survival/kaplan_meier_by_cluster.png)

*Kaplan-Meier survival curves showing distinct outcomes between proteomic-defined tumor subgroups.*

Proteomic and transcriptomic integration reveals biologically distinct tumor subgroups in clear cell renal cell carcinoma

📌 Overview

This project explores clear cell renal cell carcinoma (ccRCC) using a multi-omics approach integrating:

Clinical survival data
RNA-seq (transcriptomics)
RPPA (proteomics)

The goal was not to build a complex pipeline for its own sake, but to answer a clinically relevant question:

Can molecular patterns derived from proteomic data identify patient subgroups with different survival outcomes, and are these patterns supported at the transcriptomic level?

🧠 Clinical Context

ccRCC is a biologically heterogeneous disease.
Patients with similar clinical staging often show markedly different outcomes.

Understanding molecular subtypes linked to prognosis is essential for:

Risk stratification
Therapeutic decision-making
Precision oncology

⚙️ Data Sources

Data were obtained from the TCGA PanCancer Atlas, including:

RNA-seq expression data
RPPA proteomic profiles
Curated clinical survival outcomes (OS, OS.time)

Final cohort:

537 patients (clinical)
475 patients with complete multi-omics data

🔬 Methodology

1. Data Preparation
Extraction of KIRC samples from PanCancer dataset
Harmonization of patient identifiers across datasets
Construction of matched multi-omics cohort
2. Proteomic Analysis (RPPA)
Patient-level aggregation
Missing value imputation (median)
Variance filtering
PCA for dimensionality reduction
3. Patient Clustering
K-means clustering (k = 2) on PCA components
Identification of two proteomic subgroups
4. Survival Analysis
Kaplan–Meier survival curves
Comparison of overall survival between clusters
5. Differential Expression Analysis
Wilcoxon rank-sum test (RPPA and RNA)
FDR correction
Identification of top features per cluster
6. Multi-Omics Integration
Mapping proteins to gene symbols
Cross-layer concordance analysis
Identification of markers consistently deregulated at both protein and RNA levels

📊 Key Results

🔹 Patient Stratification

Two proteomic clusters were identified:

Cluster 1: 86 patients
Cluster 2: 389 patients
🔹 Survival Differences

Cluster 1 showed a higher event rate, suggesting worse prognosis:

Cluster 1 → ~39% events
Cluster 2 → ~34% events

Although moderate, the difference is consistent and biologically meaningful.

🔹 Multi-Omics Concordance

A set of 17 concordant markers was identified across proteomic and transcriptomic layers.

Key pathways enriched in the high-risk group:

🧬 Cell Cycle Activation
CCNE1 (Cyclin E1)
CCNB1 (Cyclin B1)
🧬 DNA Damage Response
CHEK2
RAD51 (proteomic layer)
🧬 Apoptosis Regulation
BAK1
BCL2L1 (BCL-XL)
DIABLO (SMAC)
🧬 Growth & Metabolism
IGFBP2
ASNS
TFRC

🧠 Biological Interpretation

The high-risk cluster is characterized by:

Increased proliferative signaling
Enhanced DNA repair activity
Dysregulated apoptosis
Metabolic adaptation

This profile is consistent with a more aggressive tumor phenotype, providing a mechanistic explanation for the observed survival differences.

🧾 Key Takeaway

Multi-omics integration of TCGA-KIRC reveals a biologically distinct tumor subgroup characterized by coordinated activation of cell cycle, DNA damage response, and apoptosis pathways, associated with worse survival outcomes.

📁 Project Structure

scripts/
  01_download_tcga_data.R
  02_clean_clinical_survival.R
  03_build_rnaseq_matrix.R
  04_match_multiomics_cohort.R
  05_rppa_pca_clustering.R
  08_survival_by_cluster.R
  09_rppa_differential_analysis.R
  10_rnaseq_differential_by_cluster.R
  11_multiomics_integration.R

data/
  raw/
  processed/

results/
  clustering/
  survival/
  differential/
  integration/

🧪 Reproducibility

All analyses were performed using:

R (4.5.x)
Packages:
dplyr
readr
ggplot2
survival

Each step is modular and reproducible via standalone scripts.

Project Value

This project demonstrates:

Real multi-omics integration (not simulated)
Clinically meaningful survival analysis
Cross-layer biological validation
Structured and reproducible workflow

👨‍⚕️ Author

Cristian Arias, MD
Nephrologist | Healthcare Data Scientist | Bioinformatics 


Future Directions

Pathway enrichment analysis (GO / KEGG)
External validation cohorts
Integration with mutational data
Development of prognostic risk scores