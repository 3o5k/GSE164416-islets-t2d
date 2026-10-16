# 02_qc.R
# Low-count filtering, variance-stabilising transformation, PCA, and an
# exocrine contamination score that is added to the DESeq2 design.
# Input : results/counts_sub.rds, results/coldata.rds
# Output: results/dds_filtered.rds, results/sample_qc.csv,
#         figures/pca_condition_batch.png, figures/pca_exocrine.png


library(DESeq2)
library(tidyverse)

counts_sub <- readRDS("results/counts_sub.rds")
coldata    <- readRDS("results/coldata.rds")

# Recode batch dates to safe, chronologically sortable labels (e.g. 2014_01_16)
coldata$batch <- factor(format(as.Date(as.character(coldata$batch), format = "%d/%m/%Y"),
                               "%Y_%m_%d"))
levels(coldata$batch)

# 1. Build the DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                              colData   = coldata,
                              design    = ~ batch + condition)
dds

# 2. Keep genes with >= 10 counts in at least as many samples as the smallest group
keep_genes <- rowSums(counts(dds) >= 10) >= min(table(dds$condition))
table(keep_genes)
dds <- dds[keep_genes, ]
dim(dds)

# 3. Variance-stabilising transformation (for QC and plots only)
vsd <- vst(dds, blind = TRUE)

# 4. PCA on the 500 most variable genes
pca_data    <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = condition, shape = batch)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance"),
       title = "PCA of islet transcriptomes (VST)") +
  theme_bw()

p_pca
ggsave("figures/pca_condition_batch.png", p_pca, width = 7, height = 5, dpi = 300)

# 5. What drives PC1? Correlate PCs with endocrine and exocrine markers
markers <- c(INS = "3630", GCG = "2641",
             PRSS1 = "5644", CPA1 = "1357", CELA3A = "10136")
marker_expr <- t(assay(vsd)[markers, ])
colnames(marker_expr) <- names(markers)
round(cor(pca_data[, c("PC1", "PC2")], marker_expr), 2)

# 6. Identify the sample far to the right on PC1
pca_data[pca_data$PC1 > 50, ]

# 7. Exocrine contamination score: mean VST of acinar marker genes
exo_ids <- c(PRSS1 = "5644", CPA1 = "1357", CELA3A = "10136",
             CTRB1 = "1504", CPB1 = "1360", PNLIP = "5406")
exo_ids <- exo_ids[exo_ids %in% rownames(vsd)]

qc <- as.data.frame(colData(dds))
qc$exocrine <- colMeans(assay(vsd)[exo_ids, ])

# Does contamination differ between ND and T2D?
tapply(qc$exocrine, qc$condition, median)
wilcox.test(exocrine ~ condition, data = qc)

# The five most contaminated samples
head(qc[order(qc$exocrine, decreasing = TRUE), c("condition", "batch", "exocrine")], 5)

# PCA coloured by the exocrine score
pca_data$exocrine <- qc$exocrine

p_pca_exo <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = exocrine, shape = condition)) +
  geom_point(size = 3) +
  scale_colour_viridis_c() +
  labs(x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance"),
       colour = "Exocrine\nscore",
       title = "PC1 tracks exocrine contamination") +
  theme_bw()

p_pca_exo
ggsave("figures/pca_exocrine.png", p_pca_exo, width = 7, height = 5, dpi = 300)

# 8. Add the scaled exocrine score as a covariate in the design
dds$exocrine <- as.numeric(scale(qc$exocrine))
design(dds) <- ~ batch + exocrine + condition

# 9. Save QC outputs
write.csv(qc, "results/sample_qc.csv")
saveRDS(dds, "results/dds_filtered.rds")