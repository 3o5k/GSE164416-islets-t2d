# 04_gsea.R
# Pre-ranked gene set enrichment analysis (fgsea) with MSigDB Hallmark sets,
# genes ranked by the DESeq2 Wald statistic.
# Input : results/deseq2_T2D_vs_ND.csv
# Output: results/gsea_hallmark_T2D_vs_ND.csv,
#         figures/gsea_hallmark.png, figures/gsea_oxphos.png

library(tidyverse)
library(fgsea)
library(msigdbr)

res_df <- read.csv("results/deseq2_T2D_vs_ND.csv")

# 1. Ranked gene list by the Wald statistic (one value per gene symbol)
ranks <- res_df %>%
  filter(!is.na(stat), !is.na(Symbol)) %>%
  group_by(Symbol) %>%
  summarise(stat = stat[which.max(abs(stat))]) %>%
  deframe()
ranks <- sort(ranks, decreasing = TRUE)
length(ranks)

# 2. MSigDB Hallmark gene sets
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
length(pathways)

# 3. Run GSEA (nproc = 1 keeps it single-threaded and reproducible)
set.seed(42)
gsea_res <- fgsea(pathways = pathways, stats = ranks,
                  minSize = 15, maxSize = 500, nproc = 1) %>%
  arrange(padj)

head(gsea_res %>% select(pathway, NES, padj, size), 20)

gsea_res %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
  write.csv("results/gsea_hallmark_T2D_vs_ND.csv", row.names = FALSE)

# 4. Bar plot of significant Hallmark pathways
gsea_plot <- gsea_res %>%
  filter(padj < 0.05) %>%
  mutate(pathway   = str_remove(pathway, "HALLMARK_") %>% str_replace_all("_", " "),
         direction = ifelse(NES > 0, "Up in T2D", "Down in T2D"))
nrow(gsea_plot)

p_gsea <- ggplot(gsea_plot, aes(x = NES, y = reorder(pathway, NES), fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Down in T2D" = "#2166AC", "Up in T2D" = "#B2182B")) +
  labs(x = "Normalised enrichment score (NES)", y = NULL, fill = NULL,
       title = "Hallmark pathways, T2D vs ND islets (FDR < 0.05)") +
  theme_bw()

p_gsea
ggsave("figures/gsea_hallmark.png", p_gsea, width = 8, height = 7, dpi = 300)

# 5. Enrichment plot for oxidative phosphorylation
p_oxphos <- plotEnrichment(pathways[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks) +
  labs(title = "Oxidative phosphorylation (down in T2D)")
p_oxphos
ggsave("figures/gsea_oxphos.png", p_oxphos, width = 7, height = 4, dpi = 300)