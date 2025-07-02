library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)

# Load data
count_df <- read.csv('/data/haochun/fyp/TCGA_diff_int/ESCA_STAD_gene', sep = '\t', header = TRUE, row.names = 1)
meta_df <- read.csv('/data/haochun/fyp/TCGA_diff_int/ESCA_STAD', sep = '\t', header = TRUE, row.names = 1)

#Ensure sample names match
common_samples <- intersect(colnames(count_df), rownames(meta_df))
count_df <- count_df[, common_samples]
meta_df <- meta_df[common_samples, , drop = FALSE]

#Filter low-expressed genes 
count_df <- count_df[rowSums(count_df) >= 10, ]

# Create DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = count_df,
                              colData = meta_df,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ] 

# Run DE analysis
dds <- DESeq(dds)

#  Get results with defined contrast 
res <- results(dds, contrast = c("condition", "STAD", "ESCA"))
res <- lfcShrink(dds, coef = "condition_STAD_vs_ESCA", type = "apeglm")  # 更稳健的 log2FC

# Clean results 
res_df <- as.data.frame(res)
res_df <- res_df %>% drop_na()

# Define thresholds 
padj_cutoff <- 0.05
log2fc_cutoff <- log2(2)

# Filter DEGs 
deg_df <- res_df %>%
  filter(padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)

# Save results 
write.csv(res_df, "/data/haochun/fyp/TCGA_diff_int/ESCA_STAD_all_DESeq2.csv")
write.csv(deg_df, "/data/haochun/fyp/TCGA_diff_int/ESCA_STAD_diffsig.csv")

#  QC: PCA & Heatmap 
rld <- rlog(dds, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of samples")

top_genes <- head(order(res_df$padj), 50)
pheatmap(assay(rld)[top_genes, ], 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_col = meta_df, show_rownames = FALSE)

upregulated <- deg_df %>% filter(log2FoldChange > 0)
downregulated <- deg_df %>% filter(log2FoldChange < 0)
write.csv(upregulated, "/data/haochun/fyp/TCGA_diff_int/ESCA_STAD_up.csv")
write.csv(downregulated, "/data/haochun/fyp/TCGA_diff_int/ESCA_STAD_down.csv")
