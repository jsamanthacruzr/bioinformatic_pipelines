# ==============================================================================
# RNA-seq Differential Expression Analysis (DESeq2)
# Target: ASOs / Treatment Project (DMSO vs TPL conditions)
# 
# References:
# - Main workflow: Sboner et al., Clinical Genomics Lecture
#   (https://physiology.med.cornell.edu/faculty/mason/lab/clinicalgenomics/lectures/2015/Lecture20_Sboner.pdf)
# - Bioconductor DESeq2 Vignette:
#   (https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
# - Computational Biology R Presentation:
#   (https://lashlock.github.io/compbio/R_presentation.html)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(clusterProfiler) # For GO and KEGG enrichment analysis
library(edgeR)
library(biomaRt)
library(EnhancedVolcano)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)

# ------------------------------------------------------------------------------
# 2. Environment Setup & Relative Paths
# ------------------------------------------------------------------------------
# Set working directory to project root (uncomment and adjust if needed)
# setwd("./")

results_dir <- "./results/"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 3. Load & Filter Count Data
# ------------------------------------------------------------------------------
# Load raw read count table
countTable_deseq <- read.csv("star_counts.csv", header = TRUE, row.names = 1)
head(countTable_deseq)

# Define experimental conditions (2 groups, 3 replicates each)
conditions <- c("DMSO", "DMSO", "DMSO", "TPL", "TPL", "TPL")

dim(countTable_deseq)
cpmsTable <- cpm(countTable_deseq)
head(cpmsTable)

# Filter low-expressed genes using edgeR's filterByExpr
keep.exprs <- filterByExpr(countTable_deseq, group = conditions)
countTableFilter_deseq <- countTable_deseq[keep.exprs, ]
dim(countTableFilter_deseq)
head(countTableFilter_deseq)

# ------------------------------------------------------------------------------
# 4. Experimental Metadata & DESeq2 Matrix Setup
# ------------------------------------------------------------------------------
# Create sample metadata data frame
colData <- DataFrame(condition = factor(c("DMSO", "DMSO", "DMSO", "TPL", "TPL", "TPL")))

# Set 'DMSO' as the baseline reference level for comparison (DMSO vs TPL)
colData$condition <- relevel(colData$condition, ref = "DMSO")

# Create DESeq2 dataset from count matrix
dds <- DESeqDataSetFromMatrix(
  countData = countTableFilter_deseq,
  colData   = colData,
  design    = ~ condition
)

# ------------------------------------------------------------------------------
# 5. Run Differential Expression & Quality Control (PCA)
# ------------------------------------------------------------------------------
# Run DESeq pipeline
dds <- DESeq(dds)

# Regularized log transformation for visualization
rld <- rlog(dds)

# Principal Component Analysis (PCA) plot
plotPCA(rld)

# ------------------------------------------------------------------------------
# 6. Extract & Annotate Results (biomaRt)
# ------------------------------------------------------------------------------
# Extract DESeq2 results table
res_deseq <- results(dds)
head(results(dds))
summary(res_deseq)

# Sort results by log2 fold change
res_deseqOrdered <- res_deseq[order(res_deseq$log2FoldChange), ]
head(res_deseqOrdered)

# Store Ensembl IDs from rownames
res_deseqOrdered$ensemblid <- rownames(res_deseqOrdered)

# Annotate with HGNC symbols via biomaRt
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
res_deseq_gene <- getBM(
  attributes = c('hgnc_symbol', 'ensembl_gene_id'),
  filters    = c('ensembl_gene_id'),
  values     = gsub("\\..*", "", row.names(res_deseqOrdered)),
  mart       = mart
)

# Clean Ensembl version suffixes (e.g., ENSG00000123456.1 -> ENSG00000123456)
res_deseqOrdered$ensemblid <- gsub("\\..*", "", res_deseqOrdered$ensemblid)

# Merge gene symbols with DESeq2 results
res_deseqOrdered <- merge(
  x      = res_deseq_gene, 
  by.x   = "ensembl_gene_id", 
  y      = as.data.frame(res_deseqOrdered), 
  by.y   = "ensemblid", 
  all    = TRUE
)

# ------------------------------------------------------------------------------
# 7. Filter Significant Differentially Expressed Genes (DEGs)
# Criteria: adjusted p-value < 0.05 & |log2FoldChange| >= 1.5
# ------------------------------------------------------------------------------
sig <- res_deseqOrdered[
  !is.na(res_deseqOrdered$padj) & 
    res_deseqOrdered$padj < 0.05 & 
    abs(res_deseqOrdered$log2FoldChange) >= 1.5, 
]

head(sig)
dim(sig)

# Separate downregulated and upregulated genes
keepDown_deseq <- sig$log2FoldChange < 0
dim(sig[keepDown_deseq, ])

keepUp_deseq <- sig$log2FoldChange > 0
dim(sig[keepUp_deseq, ])

# Add regulation status column
sig$expression <- ""
sig[keepDown_deseq, "expression"] <- "Down"
sig[keepUp_deseq, "expression"]   <- "Up"
head