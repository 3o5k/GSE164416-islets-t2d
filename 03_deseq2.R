# 03_deseq2.R
# DESeq2 model (~ batch + exocrine + condition), apeglm LFC shrinkage,
# gene annotation and volcano plot for T2D vs ND.
# Input : results/dds_filtered.rds, data/Human.GRCh38.p13.annot.tsv.gz
# Output: results/dds_deseq.rds, results/deseq2_T2D_vs_ND.csv,
#         figures/dispersion_estimates.png, figures/volcano_T2D_vs_ND.png

library(DESeq2)
library(tidyverse)
library(ggrepel)

dds <- readRDS("results/dds_filtered.rds")

# 1. Fit the model: size factors, dispersions, negative binomial GLM, Wald test
dds <- DESeq(dds)
saveRDS(dds, "results/dds_deseq.rds")
resultsNames(dds)

# 2. Results: T2D vs ND
res <- results(dds, name = "condition_T2D_vs_ND", alpha = 0.05)
summary(res)

# 3. Model diagnostic: dispersion estimates
plotDispEsts(dds)

png("figures/dispersion_estimates.png", width = 7, height = 5, units = "in", res = 300)
plotDispEsts(dds)
dev.off()

# 4. Shrink log2 fold changes (apeglm) for ranking and plotting
res_shr <- lfcShrink(dds, coef = "condition_T2D_vs_ND", type = "apeglm")

# 5. Add gene symbols from the NCBI annotation file
annot <- read.delim("data/Human.GRCh38.p13.annot.tsv.gz", check.names = FALSE)
annot <- annot %>%
  mutate(GeneID = as.character(GeneID)) %>%
  select(GeneID, Symbol, Description)

res_df <- as.data.frame(res) %>%
  rownames_to_column("GeneID") %>%
  mutate(log2FC_shrunk = res_shr$log2FoldChange) %>%
  left_join(annot, by = "GeneID") %>%
  arrange(padj)

# 6. Top genes and a stricter count (FDR < 0.05 and |shrunken LFC| > 1)
head(res_df %>% select(Symbol, baseMean, log2FC_shrunk, padj), 20)
sum(res_df$padj < 0.05 & abs(res_df$log2FC_shrunk) > 1, na.rm = TRUE)

write.csv(res_df, "results/deseq2_T2D_vs_ND.csv", row.names = FALSE)

# 7. Volcano plot
volc <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(status = case_when(
    padj < 0.05 & log2FC_shrunk >  1 ~ "Up in T2D",
    padj < 0.05 & log2FC_shrunk < -1 ~ "Down in T2D",
    TRUE                             ~ "Not significant"))

table(volc$status)

top_labels <- volc %>%
  filter(status != "Not significant") %>%
  slice_min(padj, n = 15)

p_volc <- ggplot(volc, aes(x = log2FC_shrunk, y = -log10(padj), colour = status)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_text_repel(data = top_labels, aes(label = Symbol),
                  size = 3, colour = "black", max.overlaps = 20) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("Down in T2D"     = "#2166AC",
                                 "Not significant" = "grey75",
                                 "Up in T2D"       = "#B2182B")) +
  labs(x = "Shrunken log2 fold change (T2D vs ND)",
       y = "-log10 adjusted p-value",
       title = "Differential expression in T2D islets",
       colour = NULL) +
  theme_bw()

p_volc
ggsave("figures/volcano_T2D_vs_ND.png", p_volc, width = 7, height = 6, dpi = 300)

# 8. Are the up-regulated genes mostly low-expressed?
volc %>%
  group_by(status) %>%
  summarise(n = n(), median_baseMean = median(baseMean))