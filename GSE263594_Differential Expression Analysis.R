# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# Set working directory (modify as needed)
setwd("C:/Users/User/Downloads/GSE263594_raw_count_matrix")  # Update this path

# Load raw count matrix
counts <- read.csv("GSE263594_raw_count_matrix.csv", header=TRUE, row.names=1)

# Load metadata
metadata <- read.csv("meta_data.csv", header=TRUE, row.names=1)

# Convert condition to a factor
metadata$Condition <- factor(metadata$Condition, levels = c("Control", "Treated"))

dim(counts)  # Check the dimensions (genes x samples)
dim(metadata)  # Ensure the sample names match count matrix columns

all(rownames(metadata) %in% colnames(counts))  # Should return TRUE
colnames(counts) <- gsub("\\.csv", "", colnames(counts))  # Remove file extensions if needed
rownames(metadata) <- gsub("\\.csv", "", rownames(metadata))  # Remove file extensions
counts <- round(counts)  # Convert floating-point values to nearest integers

summary(counts)  # Should now show only whole numbers
metadata$Condition <- factor(metadata$Condition, levels = c("Control", "Treated", "Alternate Control"))
levels(metadata$Condition)
str(metadata)
all(rownames(metadata) %in% colnames(counts))  # Should return TRUE


# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Condition)

# Run DESeq2 normalization & differential expression analysis
dds <- DESeq(dds)

res <- results(dds, contrast = c("Condition", "Treated", "Control"))


# Convert results to a dataframe
res_df <- as.data.frame(res)

# Add Gene IDs as a separate column
res_df$Gene_ID <- rownames(res_df)

# Select relevant columns
res_df <- res_df[, c("Gene_ID", "log2FoldChange", "pvalue", "padj")]

# Classify genes as Upregulated, Downregulated, or Not Significant
res_df$Regulation <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Upregulated",
                            ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Downregulated", "Not Significant"))

# Save results to a CSV file
write.csv(res_df, "DGE_results_GSE263594.csv", row.names = FALSE)

# Print summary of differential expression
summary(res)

library(EnhancedVolcano)

EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "pvalue",
                title = "Differential Gene Expression (GSE263594)",
                pCutoff = 0.05,
                FCcutoff = 1)

library(pheatmap)

# Select top 50 differentially expressed genes
top_genes <- rownames(res)[order(res$padj, na.last = NA)][1:50]

# Generate heatmap
pheatmap(assay(vst(dds))[top_genes,], cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         main = "Heatmap of Top 50 Differentially Expressed Genes")

# Classify genes based on log2FoldChange and adjusted p-value (FDR)
res_df$Regulation <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Upregulated",
                            ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Downregulated", "Not Significant"))

# View first few rows
head(res_df)

upregulated_genes <- res_df[res_df$Regulation == "Upregulated", ]
write.csv(upregulated_genes, "Upregulated_Genes.csv", row.names = FALSE)

downregulated_genes <- res_df[res_df$Regulation == "Downregulated", ]
write.csv(downregulated_genes, "Downregulated_Genes.csv", row.names = FALSE)
table(res_df$Regulation)

# Top 10 Upregulated
head(upregulated_genes[order(-upregulated_genes$log2FoldChange), ], 10)

# Top 10 Downregulated
head(downregulated_genes[order(downregulated_genes$log2FoldChange), ], 10)

library(EnhancedVolcano)

EnhancedVolcano(res_df,
                lab = res_df$Gene_ID,  
                x = "log2FoldChange",
                y = "pvalue",
                title = "Volcano Plot of Differential Expression",
                pCutoff = 0.05,
                FCcutoff = 1,
                labSize = 3.0)

head(res_df)
summary(res_df)



#extra part
res <- results(dds, contrast = c("Condition", "Treated", "Control"), independentFiltering = FALSE)
summary(res)

res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$padj), ]
res_df$Regulation <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 0, "Upregulated",
                            ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < 0, "Downregulated", "Not Significant"))

write.csv(res_df, "DEG_results_final.csv", row.names = TRUE)
head(res_df[order(-res_df$log2FoldChange), ], 10)
head(res_df[order(res_df$log2FoldChange), ], 10)

library(ggplot2)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value")








