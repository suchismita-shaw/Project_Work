**RNA-Seq Differential Expression Analysis using DESeq2**

📌 **Overview**
This repository contains R scripts and datasets for performing differential expression analysis of RNA-seq data using DESeq2. The analysis identifies upregulated and downregulated genes based on given conditions.
 Data
 
**Raw Count Matrix**: GSE263594_raw_count_matrix.csv

**Metadata (Sample Information)**: meta_data.csv

Differential Expression Results: DEG_results_final.csv, Significant_DEGs.csv
nstallation & Dependencies
1️⃣ Install Required R Packages
Make sure you have R and the required packages installed:
install.packages(c("DESeq2", "ggplot2", "pheatmap", "dplyr"))
Alternatively, use:
install.packages("BiocManager")
BiocManager::install("DESeq2")
The analysis follows these steps:

1️⃣ Preprocessing the Data

Load raw count matrix (GSE263594_raw_count_matrix.csv)
Load sample metadata (meta_data.csv)
Ensure row and column names match
2️⃣ Differential Expression Analysis with DESeq2

Create DESeq2 dataset
Normalize and perform differential expression analysis
Extract upregulated and downregulated genes
3️⃣ Results & Visualization

Save results (DEG_results_final.csv)
Generate a Volcano Plot and MA Plot for visualization
Example Results
head(res_df[order(-res_df$log2FoldChange), ], 10)
Top Downregulated Genes
head(res_df[order(res_df$log2FoldChange), ], 10)
Volcano Plot
library(ggplot2)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value")
Contact
If you have any questions or suggestions, feel free to open an issue or reach out!

📧 Email: suchismita10.shaw@gmail.com
