# ==============================================================================
# Multi-omics Cross-Comparison Analysis Pipeline
# Comparative Scatter Plots:
#   1. PRO-seq Pause Region vs. PRO-seq Gene Body
#   2. RNA-seq vs. PRO-seq (Gene Body)
#   3. CUT&Tag (H3K27ac) vs. PRO-seq (Pause Region)
#   4. CUT&Tag (H3K27ac) vs. RNA-seq
#
# Author: Samantha Cruz (2024)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(edgeR)
library(biomaRt)
library(EnhancedVolcano)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(dplyr)

# ------------------------------------------------------------------------------
# 2. Environment Setup & Relative Paths
# ------------------------------------------------------------------------------
# Set working directory to project root (uncomment and adjust if needed)
# setwd("./")

results_dir <- "./results/"
data_dir    <- "./data/"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# ==============================================================================
# SECTION 1: PRO-seq Pause Region vs. PRO-seq Gene Body
# ==============================================================================

# Load pause region dataset and collapse isoforms by gene name
pause <- read.csv(file.path(data_dir, "pause_counts_TPL1uM_vs_DMSO_All.tsv"), header = TRUE, sep = "\t")
pause$GENENAME <- gsub('"', '', trimws(pause$GENENAME))
pause <- pause %>%
  group_by(GENENAME) %>%
  slice(1) %>%
  ungroup()

# Load gene body dataset and clean gene names
gene_body <- read.csv(file.path(data_dir, "counts_TPL1uM_vs_DMSO_All.tsv"), header = TRUE, sep = "\t")
gene_body$GENENAME <- gsub('"', '', trimws(gene_body$GENENAME))

# Merge datasets and format column names
data <- merge(x = pause, by.x = "GENENAME", y = gene_body, by.y = "GENENAME", all = TRUE)

data_clean <- data[, c(1, 8, 12, 20, 24)] 
names(data_clean)[names(data_clean) == "GENENAME"]       <- "gene_names"
names(data_clean)[names(data_clean) == "log2FoldChange.x"] <- "pause_FC"
names(data_clean)[names(data_clean) == "padj.x"]           <- "p_values_pause"
names(data_clean)[names(data_clean) == "log2FoldChange.y"] <- "gene_body_FC"
names(data_clean)[names(data_clean) == "padj.y"]           <- "p_values_gene_body"

# Thresholds
fold_change_threshold <- 1.0
p_value_threshold     <- 0.05

# Categorize genes by expression pattern
data_clean$Color <- "grey"

# 1. Red: Upregulated in both Pause and Gene Body
data_clean$Color[data_clean$pause_FC > fold_change_threshold & 
                   data_clean$gene_body_FC > fold_change_threshold & 
                   data_clean$p_values_pause < p_value_threshold & 
                   data_clean$p_values_gene_body < p_value_threshold] <- "#DC3440"

# 2. Blue: Downregulated in both Pause and Gene Body
data_clean$Color[data_clean$pause_FC < -fold_change_threshold & 
                   data_clean$gene_body_FC < -fold_change_threshold & 
                   data_clean$p_values_pause < p_value_threshold & 
                   data_clean$p_values_gene_body < p_value_threshold] <- "#3690AE"

# 3. Green: Downregulated in Pause, Upregulated in Gene Body
data_clean$Color[data_clean$pause_FC < -fold_change_threshold & 
                   data_clean$gene_body_FC > fold_change_threshold & 
                   data_clean$p_values_pause < p_value_threshold & 
                   data_clean$p_values_gene_body < p_value_threshold] <- "#55A154"

# 4. Orange: Upregulated in Pause, Downregulated in Gene Body
data_clean$Color[data_clean$pause_FC > fold_change_threshold & 
                   data_clean$gene_body_FC < -fold_change_threshold & 
                   data_clean$p_values_pause < p_value_threshold & 
                   data_clean$p_values_gene_body < p_value_threshold] <- "#F28227"

# Scatter plot
plot(data_clean$pause_FC, data_clean$gene_body_FC, 
     xlab = "Log2 Fold Change in Pause Regions", 
     ylab = "Log2 Fold Change in Gene Body Regions", 
     main = "Pause vs Gene Body - 4h 250nM TPL", 
     pch  = 16, col = data_clean$Color)

text(data_clean$pause_FC, data_clean$gene_body_FC, labels = data_clean$gene_names, pos = 4, cex = 0.8, col = "black")
abline(v = c(-1, 1), col = "black", lty = 2)
abline(h = c(-1, 1), col = "black", lty = 2)

# Export processed data and calculate correlation
expression_counts <- data_clean %>% count(Color)
print(expression_counts)

data_editNA <- na.omit(data_clean[, c("pause_FC", "gene_body_FC")])
write.table(data_clean, file.path(results_dir, "comparison_pause-vs-genebody-1h1uM.csv"), row.names = FALSE, col.names = TRUE, sep = ",")

correlation_test <- cor.test(data_editNA$pause_FC, data_editNA$gene_body_FC)
cat("Pause vs Gene Body Correlation:\n")
print(correlation_test$estimate)
print(correlation_test$p.value)

# ==============================================================================
# SECTION 2: RNA-seq vs. PRO-seq (Gene Body)
# ==============================================================================

# Load RNA-seq DESeq2 results
RNAseq <- read.csv(file.path(data_dir, "All_Results_DESeq2_K562_htseqcounts.csv"), header = TRUE, sep = ",")
RNAseq$ensembl_gene_id <- gsub('"', '', trimws(RNAseq$ensembl_gene_id))

# Merge RNA-seq and PRO-seq Gene Body datasets
data <- merge(x = RNAseq, by.x = "ensembl_gene_id", y = gene_body, by.y = "GENENAME", all = TRUE)

data_clean <- data[, c(1, 4, 8, 22, 26)] 
names(data_clean)[names(data_clean) == "ensembl_gene_id"] <- "gene_names"
names(data_clean)[names(data_clean) == "log2FoldChange.x"] <- "RNAseq_FC"
names(data_clean)[names(data_clean) == "padj.x"]           <- "p_values_RNAseq"
names(data_clean)[names(data_clean) == "log2FoldChange.y"] <- "gene_body_FC"
names(data_clean)[names(data_clean) == "padj.y"]           <- "p_values_gene_body"

# Thresholds
fold_change_threshold_PROseq <- 1.0
fold_change_threshold_RNAseq <- 1.5
p_value_threshold            <- 0.05

# Categorize genes
data_clean$Color <- "grey"

# 1. Red: Upregulated in both RNA-seq and PRO-seq Gene Body
data_clean$Color[data_clean$RNAseq_FC > fold_change_threshold_RNAseq & 
                   data_clean$gene_body_FC > fold_change_threshold_PROseq & 
                   data_clean$p_values_RNAseq < p_value_threshold & 
                   data_clean$p_values_gene_body < p_value_threshold] <- "#DC3440"

# 2. Blue: Downregulated in both RNA-seq and PRO-seq Gene Body
data_clean$Color[data_clean$RNAseq_FC < -fold_change_threshold_RNAseq & 
                   data_clean$gene_body_FC < -fold_change_threshold_PROseq & 
                   data_clean$p_values_RNAseq < p_value_threshold & 
                   data_clean$p_values_gene_body < p_value_threshold] <- "#3690AE"

# 3. Green: Downregulated in RNA-seq, Upregulated in PRO-seq Gene Body
data_clean$Color[data_clean$RNAseq_FC < -fold_change_threshold_RNAseq & 
                   data_clean$gene_body_FC > fold_change_threshold_PROseq & 
                   data_clean$p_values_RNAseq < p_value_threshold & 
                   data_clean$p_values_gene_body < p_value_threshold] <- "#55A154"

# 4. Orange: Upregulated in RNA-seq, Downregulated in PRO-seq Gene Body
data_clean$Color[data_clean$RNAseq_FC > fold_change_threshold_RNAseq & 
                   data_clean$gene_body_FC < -fold_change_threshold_PROseq & 
                   data_clean$p_values_RNAseq < p_value_threshold & 
                   data_clean$p_values_gene_body < p_value_threshold] <- "#F28227"

# Scatter plot
plot(data_clean$RNAseq_FC, data_clean$gene_body_FC, 
     xlab = "Log2 Fold Change in RNAseq", 
     ylab = "Log2 Fold Change in PROseq", 
     main = "RNAseq vs PROseq - 4h 250nM TPL", 
     pch  = 16, col = data_clean$Color)

text(data_clean$RNAseq_FC, data_clean$gene_body_FC, labels = data_clean$gene_names, pos = 4, cex = 0.8, col = "black")
abline(v = c(-1.5, 1.5), col = "black", lty = 2)
abline(h = c(-1.0, 1.0), col = "black", lty = 2)

# Export processed data and calculate correlation
expression_counts <- data_clean %>% count(Color)
print(expression_counts)

data_editNA <- na.omit(data_clean[, c("RNAseq_FC", "gene_body_FC")])
write.table(data_clean, file.path(results_dir, "comparison_RNAseq-vs-genebody-1h1uM.csv"), row.names = FALSE, col.names = TRUE, sep = ",")

correlation_test <- cor.test(data_editNA$RNAseq_FC, data_editNA$gene_body_FC)
cat("RNA-seq vs PRO-seq Gene Body Correlation:\n")
print(correlation_test$estimate)
print(correlation_test$p.value)

# ==============================================================================
# SECTION 3: CUT&Tag (H3K27ac) vs. PRO-seq (Pause Region)
# ==============================================================================

# Load CUT&Tag DESeq2 results
CUTTag <- read.csv(file.path(data_dir, "filtered_promoters_for-comparison_CUTTAG-K27ac.csv"), header = TRUE, sep = ",")
CUTTag$SYMBOL <- gsub('"', '', trimws(CUTTag$SYMBOL))
CUTTag <- CUTTag %>%
  group_by(SYMBOL) %>%
  slice(1) %>%
  ungroup()

# Merge CUT&Tag and PRO-seq Pause Region datasets
data <- merge(x = CUTTag, by.x = "SYMBOL", y = pause, by.y = "GENENAME", all = TRUE)

data_clean <- data[, c(1, 2, 3, 10, 14)] 
names(data_clean)[names(data_clean) == "SYMBOL"]           <- "gene_names"
names(data_clean)[names(data_clean) == "log2FoldChange.x"] <- "CUTTag_FC"
names(data_clean)[names(data_clean) == "padj.x"]           <- "p_values_CUTTag"
names(data_clean)[names(data_clean) == "log2FoldChange.y"] <- "pause_FC"
names(data_clean)[names(data_clean) == "padj.y"]           <- "p_values_pause"

# Thresholds
fold_change_threshold_PROseq <- 1.0
fold_change_threshold_CUT    <- 1.2
p_value_threshold            <- 0.05

# Categorize genes
data_clean$Color <- "grey"

# 1. Red: Upregulated in both H3K27ac and Pause Region
data_clean$Color[data_clean$CUTTag_FC > fold_change_threshold_CUT & 
                   data_clean$pause_FC > fold_change_threshold_PROseq & 
                   data_clean$p_values_CUTTag < p_value_threshold & 
                   data_clean$p_values_pause < p_value_threshold] <- "#DC3440"

# 2. Blue: Downregulated in both H3K27ac and Pause Region
data_clean$Color[data_clean$CUTTag_FC < -fold_change_threshold_CUT & 
                   data_clean$pause_FC < -fold_change_threshold_PROseq & 
                   data_clean$p_values_CUTTag < p_value_threshold & 
                   data_clean$p_values_pause < p_value_threshold] <- "#3690AE"

# 3. Green: Downregulated in H3K27ac, Upregulated in Pause Region
data_clean$Color[data_clean$CUTTag_FC < -fold_change_threshold_CUT & 
                   data_clean$pause_FC > fold_change_threshold_PROseq & 
                   data_clean$p_values_CUTTag < p_value_threshold & 
                   data_clean$p_values_pause < p_value_threshold] <- "#55A154"

# 4. Orange: Upregulated in H3K27ac, Downregulated in Pause Region
data_clean$Color[data_clean$CUTTag_FC > fold_change_threshold_CUT & 
                   data_clean$pause_FC < -fold_change_threshold_PROseq & 
                   data_clean$p_values_CUTTag < p_value_threshold & 
                   data_clean$p_values_pause < p_value_threshold] <- "#F28227"

# Scatter plot
plot(data_clean$CUTTag_FC, data_clean$pause_FC, 
     xlab = "Log2 Fold Change in H3K27ac", 
     ylab = "Log2 Fold Change in PROseq", 
     main = "H3K27ac vs PROseq - 4h 250nM TPL", 
     pch  = 16, col = data_clean$Color)

text(data_clean$CUTTag_FC, data_clean$pause_FC, labels = data_clean$gene_names, pos = 4, cex = 0.8, col = "black")
abline(v = c(-1.2, 1.2), col = "black", lty = 2)
abline(h = c(-1.0, 1.0), col = "black", lty = 2)

# Export processed data and calculate correlation
expression_counts <- data_clean %>% count(Color)
print(expression_counts)

data_editNA <- na.omit(data_clean[, c("CUTTag_FC", "pause_FC")])
write.table(data_clean, file.path(results_dir, "comparison_K27ac-vs-pause-1h1uM.csv"), row.names = FALSE, col.names = TRUE, sep = ",")

correlation_test <- cor.test(data_editNA$CUTTag_FC, data_editNA$pause_FC)
cat("H3K27ac vs PRO-seq Pause Region Correlation:\n")
print(correlation_test$estimate)
print(correlation_test$p.value)

# ==============================================================================
# SECTION 4: CUT&Tag (H3K27ac) vs. RNA-seq
# ==============================================================================

# Merge CUT&Tag and RNA-seq datasets
data <- merge(x = CUTTag, by.x = "SYMBOL", y = RNAseq, by.y = "ensembl_gene_id", all = TRUE)

data_clean <- data[, c(1, 2, 3, 6, 10)] 
names(data_clean)[names(data_clean) == "SYMBOL"]           <- "gene_names"
names(data_clean)[names(data_clean) == "log2FoldChange.x"] <- "CUTTag_FC"
names(data_clean)[names(data_clean) == "padj.x"]           <- "p_values_CUTTag"
names(data_clean)[names(data_clean) == "log2FoldChange.y"] <- "RNAseq_FC"
names(data_clean)[names(data_clean) == "padj.y"]           <- "p_values_RNAseq"

# Thresholds
fold_change_threshold_CUT    <- 1.2
fold_change_threshold_RNAseq <- 1.5
p_value_threshold            <- 0.05

# Categorize genes
data_clean$Color <- "grey"

# 1. Red: Upregulated in both H3K27ac and RNA-seq
data_clean$Color[data_clean$CUTTag_FC > fold_change_threshold_CUT & 
                   data_clean$RNAseq_FC > fold_change_threshold_RNAseq & 
                   data_clean$p_values_CUTTag < p_value_threshold & 
                   data_clean$p_values_RNAseq < p_value_threshold] <- "#DC3440"

# 2. Blue: Downregulated in both H3K27ac and RNA-seq
data_clean$Color[data_clean$CUTTag_FC < -fold_change_threshold_CUT & 
                   data_clean$RNAseq_FC < -fold_change_threshold_RNAseq & 
                   data_clean$p_values_CUTTag < p_value_threshold & 
                   data_clean$p_values_RNAseq < p_value_threshold] <- "#3690AE"

# 3. Green: Downregulated in H3K27ac, Upregulated in RNA-seq
data_clean$Color[data_clean$CUTTag_FC < -fold_change_threshold_CUT & 
                   data_clean$RNAseq_FC > fold_change_threshold_RNAseq & 
                   data_clean$p_values_CUTTag < p_value_threshold & 
                   data_clean$p_values_RNAseq < p_value_threshold] <- "#55A154"

# 4. Orange: Upregulated in H3K27ac, Downregulated in RNA-seq
data_clean$Color[data_clean$CUTTag_FC > fold_change_threshold_CUT & 
                   data_clean$RNAseq_FC < -fold_change_threshold_RNAseq & 
                   data_clean$p_values_CUTTag < p_value_threshold & 
                   data_clean$p_values_RNAseq < p_value_threshold] <- "#F28227"

# Scatter plot (Fixed Y-axis coordinate positioning for labels)
plot(data_clean$CUTTag_FC, data_clean$RNAseq_FC, 
     xlab = "Log2 Fold Change in H3K27ac", 
     ylab = "Log2 Fold Change in RNAseq", 
     main = "H3K27ac vs RNAseq", 
     pch  = 16, col = data_clean$Color)

text(data_clean$CUTTag_FC, data_clean$RNAseq_FC,