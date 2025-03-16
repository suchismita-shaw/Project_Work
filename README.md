# RNA-seq Differential Expression Analysis using DESeq2

## Overview
This repository contains the R scripts and processed data for RNA-seq differential expression analysis using DESeq2. The project aims to identify upregulated and downregulated genes between experimental conditions.

## Data Description
raw_counts.csv → Contains raw gene expression counts from RNA-seq.

metadata.csv → Contains sample conditions (e.g., Control, Treated).

DEG_results_final.csv → Contains results for all genes.

Significant_DEGs.csv → Filtered differentially expressed genes (padj < 0.05).

## Requirements
Ensure you have the following dependencies installed in R before running the analysis:
```bash
install.packages(c("ggplot2", "pheatmap", "dplyr", "DESeq2"))
```

### The analysis follows these steps:

1️⃣ Preprocessing the Data

Load raw count matrix (GSE263594_raw_count_matrix.csv)
Load sample metadata (meta_data.csv)
Ensure row and column names match

2️⃣ Differential Expression Analysis with DESeq2

Create DESeq2 dataset
Normalize and perform differential expression analysis
Extract upregulated and downregulated genes

3️⃣ Results & visualization

Save results (DEG_results_final.csv)
Generate a Volcano Plot and MA Plot for visualization

## 📊 Example Results

### 1️⃣ Top Upregulated Genes

```bash
head(res_df[order(-res_df$log2FoldChange), ], 10)
```

###  2️⃣ Top Downregulated Genes

```bash
head(res_df[order(res_df$log2FoldChange), ], 10)
```
## 📩 Contact
If you have any questions or suggestions, feel free to open an issue or reach out!

📧 Email: suchismita10.shaw@gmail.com 
