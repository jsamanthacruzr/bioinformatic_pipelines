# ==============================================================================
# CUT&Tag Data Analysis Pipeline
# Based on: https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=14
# Target: H3K27ac
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(chromVAR)        # For FRiP analysis and differential analysis
library(DESeq2)          # For differential analysis
library(ggpubr)          # For customizing figures
library(corrplot)        # For correlation matrix visualization
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v86)
library(ChIPseeker)      # For peak annotation
library(clusterProfiler) # For GO and KEGG functional enrichment
library(org.Hs.eg.db)
library(tidyverse)
library(VennDiagram)
library(ChIPseqSpikeInFree)
library(DiffBind)
library(EnhancedVolcano)

# ------------------------------------------------------------------------------
# 2. Environment Setup & Relative Paths
# ------------------------------------------------------------------------------
# Set working directory to project root (adjust relative paths as needed)
# setwd("./")

projPath <- "./default_bowtie2/"
resultsDir <- paste0(projPath, "results_R/")

if (!dir.exists(resultsDir)) {
  dir.create(resultsDir, recursive = TRUE)
}

# Define sample and histone lists
sampleList <- c("H3K27ac-DMSO_R1", "H3K27ac-DMSO_R2", "H3K27ac-DMSO_R3", 
                "H3K27ac-TPL_R1",  "H3K27ac-TPL_R2",  "H3K27ac-TPL_R3")

histList <- c("H3K27ac-DMSO", "H3K27ac-TPL")

# ------------------------------------------------------------------------------
# 3. Alignment Summary & Metrics (Bowtie2)
# ------------------------------------------------------------------------------
alignResult <- c()

for (hist in sampleList) {
  alignRes <- read.table(paste0(projPath, "bowtie2/", hist, "_bowtie2summary.txt"), header = FALSE, fill = TRUE)
  alignRate <- substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6])) - 1)
  histInfo <- strsplit(hist, "_")[[1]]
  
  alignResult <- data.frame(
    Histone            = histInfo[1], 
    Replicate          = histInfo[2], 
    SequencingDepth    = as.numeric(as.character(alignRes$V1[1])), 
    MappedFragNum_hg19 = as.numeric(as.character(alignRes$V1[4])) + as.numeric(as.character(alignRes$V1[5])), 
    AlignmentRate_hg19 = as.numeric(alignRate)
  ) %>% rbind(alignResult, .)
}

alignResult$Histone <- factor(alignResult$Histone, levels = histList)

# Export alignment summary table
write.table(alignResult, file = paste0(resultsDir, "alignResult_H3K27ac.csv"), sep = ",", row.names = FALSE)

# --- Alignment Quality Control Figures ---
fig3A <- alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth / 1e6, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "F", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Sequencing Depth (Millions)") +
  xlab("") + 
  ggtitle("A. Sequencing Depth")

fig3B <- alignResult %>% ggplot(aes(x = Histone, y = MappedFragNum_hg19 / 1e6, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "F", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Mapped Fragments (Millions)") +
  xlab("") +
  ggtitle("B. Alignable Fragments (hg19)")

fig3C <- alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate_hg19, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "F", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate (hg19)")

ggarrange(fig3A, fig3B, fig3C, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# ------------------------------------------------------------------------------
# 4. Duplication Metrics (Picard)
# ------------------------------------------------------------------------------
dupResult <- c()

for (hist in sampleList) {
  dupRes <- read.table(paste0(projPath, "picard/", hist, "_markeddupmetrics.txt"), header = TRUE, fill = TRUE)
  histInfo <- strsplit(hist, "_")[[1]]
  
  dupResult <- data.frame(
    Histone              = histInfo[1], 
    Replicate            = histInfo[2], 
    MappedFragNum_hg19   = as.numeric(as.character(dupRes$READ_PAIRS_EXAMINED[1])), 
    DuplicationRate      = as.numeric(as.character(dupRes$PERCENT_DUPLICATION[1])) * 100, 
    EstimatedLibrarySize = as.numeric(as.character(dupRes$ESTIMATED_LIBRARY_SIZE[1]))
  ) %>% 
    mutate(UniqueFragNum = MappedFragNum_hg19 * (1 - DuplicationRate / 100)) %>% 
    rbind(dupResult, .)
}

dupResult$Histone <- factor(dupResult$Histone, levels = histList)

write.table(dupResult, file = paste0(resultsDir, "dupResult_H3K27ac.csv"), sep = ",", row.names = FALSE)

# --- Duplication Plots ---
fig4A <- dupResult %>% ggplot(aes(x = Histone, y = DuplicationRate, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "F", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Duplication Rate (%)") +
  xlab("") 

fig4B <- dupResult %>% ggplot(aes(x = Histone, y = EstimatedLibrarySize, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "F", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Estimated Library Size") +
  xlab("") 

fig4C <- dupResult %>% ggplot(aes(x = Histone, y = UniqueFragNum, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "F", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Unique Fragments") +
  xlab("")

ggarrange(fig4A, fig4B, fig4C, ncol = 3, common.legend = TRUE, legend = "bottom")

# ------------------------------------------------------------------------------
# 5. Fragment Size Distribution
# ------------------------------------------------------------------------------
fragLen <- c()

for (hist in sampleList) {
  histInfo <- strsplit(hist, "_")[[1]]
  fragLen <- read.table(paste0(projPath, "fragmentLen/", hist, "_fragmentLennoduplicates.txt"), header = FALSE) %>% 
    mutate(
      fragLen    = as.numeric(V1), 
      fragCount  = as.numeric(V2), 
      Weight     = as.numeric(V2) / sum(as.numeric(V2)), 
      Histone    = histInfo[1], 
      Replicate  = histInfo[2], 
      sampleInfo = hist
    ) %>% 
    rbind(fragLen, .) 
}

fragLen$sampleInfo <- factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone    <- factor(fragLen$Histone, levels = histList)

# --- Fragment Size Figures ---
fig5A <- fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "F", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length (bp)") +
  xlab("")

fig5B <- fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Replicate)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "F") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length (bp)") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

ggarrange(fig5A, fig5B, ncol = 2)

# ------------------------------------------------------------------------------
# 6. Replicate Reproducibility Analysis
# ------------------------------------------------------------------------------
fragCount <- NULL

for (hist in sampleList) {
  if (is.null(fragCount)) {
    fragCount <- read.table(paste0(projPath, "binLen/", hist, "_bin.bed"), header = FALSE) 
    colnames(fragCount) <- c("chrom", "bin", hist)
  } else {
    fragCountTmp <- read.table(paste0(projPath, "binLen/", hist, "_bin.bed"), header = FALSE)
    colnames(fragCountTmp) <- c("chrom", "bin", hist)
    fragCount <- full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
  }
}

# Calculate Log2 Correlation Matrix
M <- cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs")  

corrplot(
  M, 
  method = "color", 
  outline = TRUE, 
  addgrid.col = "darkgray", 
  order = "hclust", 
  addrect = 3, 
  rect.col = "black", 
  rect.lwd = 3,
  cl.pos = "b", 
  tl.col = "indianred4", 
  tl.cex = 1, 
  cl.cex = 1, 
  addCoef.col = "black", 
  number.digits = 2, 
  number.cex = 1, 
  col = colorRampPalette(c("midnightblue", "white", "darkred"))(100)
)

# ------------------------------------------------------------------------------
# 7. Spike-in Free Normalization (ChIPseqSpikeInFree)
# Note: Used because no external spike-in (e.g., E. coli) was present.
# ------------------------------------------------------------------------------
metaFile <- paste0(projPath, "metadata_spikefree.txt")

bams <- c(
  "H3K27ac-DMSO_R1_noduplicatesbowtie2.bam", "H3K27ac-DMSO_R2_noduplicatesbowtie2.bam", "H3K27ac-DMSO_R3_noduplicatesbowtie2.bam",
  "H3K27ac-TPL_R1_noduplicatesbowtie2.bam",  "H3K27ac-TPL_R2_noduplicatesbowtie2.bam",  "H3K27ac-TPL_R3_noduplicatesbowtie2.bam"
)

ChIPseqSpikeInFree(bamFiles = bams, chromFile = "hg19", metaFile = metaFile, prefix = "H3K27ac")

# ------------------------------------------------------------------------------
# 8. SEACR Peak Evaluation & Reproducibility (With Spike-in Factors)
# ------------------------------------------------------------------------------
peakN <- c()
peakWidth <- c()
peakType <- c("top0.01")

for (hist in sampleList) {
  histInfo <- strsplit(hist, "_")[[1]]
  if (histInfo[1] != "IgG") {
    for (type in peakType) {
      peakInfo  <- read.table(paste0(projPath, "SEACR/", hist, "_spikein.stringent.bed"), header = FALSE, fill = TRUE) %>% mutate(width = abs(V3 - V2))
      peakN     <- data.frame(peakN = nrow(peakInfo), peakType = type, Histone = histInfo[1], Replicate = histInfo[2]) %>% rbind(peakN, .)
      peakWidth <- data.frame(width = peakInfo$width, peakType = type, Histone = histInfo[1], Replicate = histInfo[2]) %>% rbind(peakWidth, .)
    }
  }
}

# Calculate Peak Overlap Reproducibility Across Replicates
histL <- c("H3K27ac-DMSO", "H3K27ac-TPL")
repL  <- paste0("R", 1:3)
peakOverlap <- c()

for (type in peakType) {
  for (hist in histL) {
    overlap.gr <- GRanges()
    for (rep in repL) {
      peakInfo    <- read.table(paste0(projPath, "SEACR/", hist, "_", rep, "_spikein.stringent.bed"), header = FALSE, fill = TRUE)
      peakInfo.gr <- GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
      
      if (length(overlap.gr) > 0) {
        overlap.gr <- overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
      } else {
        overlap.gr <- peakInfo.gr
      }
    }
    peakOverlap <- data.frame(peakReprod = length(overlap.gr), Histone = hist, peakType = type) %>% rbind(peakOverlap, .)
  }
}

write.table(peakOverlap, paste0(resultsDir, "peakOverlap_H3K27ac_spikein.csv"), sep = ",", quote = FALSE)

peakReprod <- left_join(peakN, peakOverlap, by = c("Histone", "peakType")) %>% 
  mutate(peakReprodRate = peakReprod / peakN * 100)

write.table(peakReprod, paste0(resultsDir, "peakReprod_H3K27ac_spikein.csv"), sep = ",", quote = FALSE)

# ------------------------------------------------------------------------------
# 9. Fraction of Reads in Peaks (FRiP) Analysis
# ------------------------------------------------------------------------------
bamDir <- paste0(projPath, "bam")
inPeakData <- c()

histL <- c("H3K27ac-DMSO")
repL  <- paste0("R", 1:3)

for (hist in histL) {
  for (rep in repL) {
    peakRes <- read.table(paste0(projPath, "SEACR/", hist, "_", rep, "_spikein.stringent.bed"), header = FALSE, fill = TRUE)
    peak.gr <- GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile <- paste0(bamDir, "/", hist, "_", rep, "_noduplicatesbowtie2.bam")
    
    fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    inPeakN         <- sum(counts(fragment_counts)[, 1])
    inPeakData      <- rbind(inPeakData, data.frame(inPeakN = inPeakN, Histone = hist, Replicate = rep))
  }
}

frip <- left_join(inPeakData, alignResult, by = c("Histone", "Replicate")) %>% 
  mutate(frip = inPeakN / MappedFragNum_hg19 * 100)

# --- Summary Metrics Visualization ---
fig7A <- peakN %>% ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  facet_grid(~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") + xlab("")

fig7B <- peakWidth %>% ggplot(aes(x = Histone, y = width, fill = Histone)) +
  geom_violin() +
  facet_grid(Replicate ~ peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks (bp)") + xlab("")

fig7C <- peakReprod %>% ggplot(aes(x = Histone, y = peakReprodRate, fill = Histone, label = round(peakReprodRate, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0.1) +
  facet_grid(Replicate ~ peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") + xlab("")

fig7D <- frip %>% ggplot(aes(x = Histone, y = frip, fill = Histone, label = round(frip, 2))) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Fragments in Peaks (FRiP)") + xlab("")

ggarrange(fig7A, fig7B, fig7C, fig7D, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# ------------------------------------------------------------------------------
# 10. Differential Binding Analysis (DESeq2)
# ------------------------------------------------------------------------------
histL1 <- c("H3K27ac-DMSO")
histL2 <- c("H3K27ac-TPL")
repL1  <- paste0("R", 1:3)
repL2  <- paste0("R", 1:3)

mPeak <- GRanges()

# Construct Master Consensus Peakset
for (hist in histL1) {
  for (rep in repL1) {
    peakRes <- read.table(paste0(projPath, "SEACR/", hist, "_", rep, "_spikein.stringent.bed"), header = FALSE, fill = TRUE)
    mPeak   <- GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
  }
}

for (hist in histL2) {
  for (rep in repL2) {
    peakRes <- read.table(paste0(projPath, "SEACR/", hist, "_", rep, "_spikein.stringent.bed"), header = FALSE, fill = TRUE)
    mPeak   <- GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
  }
}

masterPeak <- reduce(mPeak)
write.table(masterPeak, paste0(resultsDir, "masterPeak_H3K27ac.csv"), sep = ",")

# Quantify Read Counts per Peak Matrix
countMat <- matrix(NA, length(masterPeak), length(histL1) * length(repL1))
i <- 1
for (hist in histL1) {
  for (rep in repL1) {
    bamFile         <- paste0(bamDir, "/", hist, "_", rep, "_noduplicatesbowtie2.bam")
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i]   <- counts(fragment_counts)[, 1]
    i               <- i + 1
  }
}
colnames(countMat) <- paste(rep(histL1, 2), rep(repL1), sep = "_")

countMat2 <- matrix(NA, length(masterPeak), length(histL2) * length(repL2))
i <- 1
for (hist in histL2) {
  for (rep in repL2) {
    bamFile          <- paste0(bamDir, "/", hist, "_", rep, "_noduplicatesbowtie2.bam")
    fragment_counts  <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat2[, i]   <- counts(fragment_counts)[, 1]
    i                <- i + 1
  }
}
colnames(countMat2) <- paste(rep(histL2, 3), rep(repL2), sep = "_")

# Combine Count Matrices
countMat <- cbind(countMat, countMat2)

# Load Normalization Factors from ChIPseqSpikeInFree
dat <- read.table(paste0(projPath, "H3K27ac_SF.txt"), sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
SF  <- dat$SF

selectR <- which(rowSums(countMat) > 5) # Filter low-count peaks
dataS   <- countMat[selectR, ]
condition <- factor(c("H3K27ac-DMSO", "H3K27ac-DMSO", "H3K27ac-DMSO", "H3K27ac-TPL", "H3K27ac-TPL", "H3K27ac-TPL"))

# DESeq2 Analysis
dds <- DESeqDataSetFromMatrix(countData = dataS, colData = DataFrame(condition), design = ~ condition)
dds <- estimateSizeFactors(dds)

# Apply Spike-in Normalization Factors
coldata <- colData(dds)
coldata$sizeFactor <- coldata$sizeFactor * SF
colData(dds) <- coldata

DDS <- DESeq(dds)
res <- results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")

countMatDiff <- cbind(dataS, res)
rownames(countMatDiff) <- selectR 

masterPeakMat       <- as.data.frame(masterPeak)
countMatDiff_peaks <- merge(masterPeakMat, countMatDiff, by = 'row.names', all.x = FALSE)
countMatDiff_peaks$Row.names <- as.numeric(countMatDiff_peaks$Row.names)
countMatDiff_peaks <- countMatDiff_peaks[order(countMatDiff_peaks$Row.names, decreasing = FALSE), ]

write.table(countMatDiff_peaks, paste0(resultsDir, "differential_expression_H3K27ac_peaks_spikeinnorm.csv"), sep = ",", quote = FALSE)

# Filter Significant Differentially Enriched Peaks (padj < 0.05 & |Log2FC| >= 1.2)
sig <- countMatDiff_peaks[!is.na(countMatDiff_peaks$padj) & 
                            countMatDiff_peaks$padj < 0.05 & 
                            abs(countMatDiff_peaks$log2FoldChange) >= 1.2, ]

keepDown_deseq <- sig$log2FoldChange < 0
keepUp_deseq   <- sig$log2FoldChange > 0

sig$expression <- ""
sig[keepDown_deseq, "expression"] <- "Down"
sig[keepUp_deseq, "expression"]   <- "Up"
sig <- sig[order(sig$log2FoldChange), ]

# Volcano Plot Visualization
EnhancedVolcano(countMatDiff_peaks, x = "log2FoldChange", y = "padj", lab = NA)

write.table(sig, paste0(resultsDir, "differential_expression_H3K27ac_peaks_significant_spikeinnorm.csv"), sep = ",", quote = FALSE)

# ------------------------------------------------------------------------------
# 11. Functional Enrichment & Peak Annotation (ChIPseeker & clusterProfiler)
# ------------------------------------------------------------------------------
samplefiles <- list.files(paste0(resultsDir, "bed"), full.names = TRUE)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("H3K27ac-DOWN", "H3K27ac-UP")

# Define Genome Annotation Database (hg19)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Peak Annotation Relative to Genomic Features
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb = txdb, 
                       tssRegion = c(-1000, 1000), annoDb = "org.Hs.eg.db", verbose = FALSE)

# Plot Feature Distribution
plotAnnoBar(peakAnnoList$`H3K27ac-DOWN`, title = "H3K27ac Downregulated Peaks - Feature Distribution") 
plotAnnoBar(peakAnnoList$`H3K27ac-UP`,   title = "H3K27ac Upregulated Peaks - Feature Distribution") 

# Plot Distance to TSS
plotDistToTSS(peakAnnoList$`H3K27ac-DOWN`, title = "H3K27ac Downregulated Peaks - Distance to TSS")
plotDistToTSS(peakAnnoList$`H3K27ac-UP`,   title = "H3K27ac Upregulated Peaks - Distance to TSS")

# Extract Annotation Data Frames
H3K27ac_DOWN_annot <- data.frame(peakAnnoList[["H3K27ac-DOWN"]]@anno)
H3K27ac_UP_annot   <- data.frame(peakAnnoList[["H3K27ac-UP"]]@anno)

write.table(H3K27ac_UP_annot, paste0(resultsDir, "H3K27ac_annot_UP_spikein.txt"), sep = "\t", quote = FALSE)

# --- Gene Ontology (GO) Enrichment Analysis ---
ego_H3K27ac_UP_ALL <- enrichGO(
  gene          = H3K27ac_UP_annot$ENSEMBL, 
  OrgDb         = org.Hs.eg.db, 
  keyType       = "ENSEMBL",
  ont           = "ALL", 
  pAdjustMethod = "BH", 
  qvalueCutoff  = 0.05, 
  readable      = TRUE
)

dotplot(ego_H3K27ac_UP_ALL, showCategory = 20)

ego_H3K27ac_DOWN_ALL <- enrichGO(
  gene          = H3K27ac_DOWN_annot$ENSEMBL, 
  OrgDb         = org.Hs.eg.db, 
  keyType       = "ENSEMBL",
  ont           = "ALL", 
  pAdjustMethod = "BH", 
  qvalueCutoff  = 0.05, 
  readable      = TRUE
)

dotplot(ego_H3K27ac_DOWN_ALL, showCategory = 20)

# Export GO Summaries
GO_summary_H3K27ac_UP <- data.frame(ego_H3K27ac_UP_ALL)
write.table(GO_summary_H3K27ac_UP, paste0(resultsDir, "clusterProfiler_GO-ALL-H3K27ac_UP_spikein.txt"), sep = "\t", quote = FALSE)

# ------------------------------------------------------------------------------
# 12. Alternative Analysis: DiffBind Pipeline (Optional)
# Reference: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/08_diffbind_differential_peaks.html
# ------------------------------------------------------------------------------
# Load Samplesheet
samples <- read.csv(paste0(projPath, "HOMER/merged/samplesheet_H3K27ac.csv"))

# Create DiffBind Object & Plot Correlation Heatmap
dbObj <- dba(sampleSheet = samples)
plot(dbObj)

# Count Reads in Consensus Peaks
dbObjcount <- dba.count(dbObj, bUseSummarizeOverlaps = TRUE)
plot(dbObjcount)

# Normalize Library Sizes
dbObjcount <- dba.normalize(dbObjcount)

# Principal Component Analysis (PCA)
dba.plotPCA(dbObjcount, DBA_FACTOR, label = DBA_CONDITION)

# Define Contrast (Baseline: DMSO)
dbObjcount <- dba.contrast(dbObjcount, reorderMeta = list(Condition = "DMSO"))

# Run Differential Analysis (DESeq2 underlying algorithm)
dbObjcount <- dba.analyze(dbObjcount)
plot(dbObjcount, contrast = 1)

# Report & Export Results
dbObjcount.DB <- dba.report(dbObjcount)
write.table(dbObjcount.DB, paste0(resultsDir, "H3K27ac_diffBind_ALL.csv"), sep = ",")

# Filter Significant Peaks (|Fold Change| >= 1.5 & FDR <= 0.05)
keepdbObj <- (dbObjcount.DB$FDR <= 0.05 & abs(dbObjcount.DB$Fold) >= 1.5)
genesED   <- dbObjcount.DB[keepdbObj, ]
genesED   <- genesED[order(genesED$Fold), ]

write.table(genesED, paste0(resultsDir, "H3K27ac_diffBind_Top.csv"), sep = ",", quote = FALSE)

# Venn Diagram of Differential Peaks
dba.plotVenn(dbObjcount, contrast = 1, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)