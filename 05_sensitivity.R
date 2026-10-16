# 05_sensitivity.R
# Re-run the full pipeline on alternative sample sets and compare with the
# main analysis: (A) all ND and T2D samples ignoring the purity filter,
# (B) the main set without the most contaminated sample (GSM5009304).
# Input : data/GSE164416_raw_counts_GRCh38.p13_NCBI.tsv.gz,
#         results/meta.rds, results/deseq2_T2D_vs_ND.csv
# Output: results/sensitivity_summary.csv

library(DESeq2)
library(tidyverse)

counts <- read.delim("data/GSE164416_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                     row.names = 1, check.names = FALSE)
meta <- readRDS("results/meta.rds")
main <- read.csv("results/deseq2_T2D_vs_ND.csv") %>%
  mutate(GeneID = as.character(GeneID))

exo_ids <- c("5644", "1357", "10136", "1504", "1360", "5406")

# Full pipeline (filter, exocrine score, DESeq2) on any chosen set of samples
run_deseq <- function(sample_ids) {
  cd <- data.frame(
    condition = factor(meta[sample_ids, "diabetes status:ch1"], levels = c("ND", "T2D")),
    batch     = factor(format(as.Date(meta[sample_ids, "library_prep_date:ch1"],
                                      format = "%d/%m/%Y"), "%Y_%m_%d")),
    row.names = sample_ids)
  dds <- DESeqDataSetFromMatrix(counts[, sample_ids], cd, design = ~ batch + condition)
  dds <- dds[rowSums(counts(dds) >= 10) >= min(table(cd$condition)), ]
  vsd <- vst(dds, blind = TRUE)
  dds$exocrine <- as.numeric(scale(colMeans(assay(vsd)[intersect(exo_ids, rownames(vsd)), ])))
  design(dds) <- ~ batch + exocrine + condition
  dds <- DESeq(dds)
  results(dds, name = "condition_T2D_vs_ND", alpha = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("GeneID")
}

# Compare an alternative analysis with the main one
compare_to_main <- function(res_alt, label) {
  m <- inner_join(main %>% select(GeneID, lfc_main = log2FoldChange, padj_main = padj),
                  res_alt %>% select(GeneID, lfc_alt = log2FoldChange, padj_alt = padj),
                  by = "GeneID")
  sig <- m %>% filter(padj_main < 0.05)
  tibble(analysis          = label,
         lfc_correlation   = cor(m$lfc_main, m$lfc_alt, use = "complete.obs"),
         sig_main          = nrow(sig),
         still_significant = sum(sig$padj_alt < 0.05, na.rm = TRUE),
         same_direction    = mean(sign(sig$lfc_main) == sign(sig$lfc_alt), na.rm = TRUE))
}

# A. All ND and T2D samples, ignoring the INS purity filter
all_ids <- rownames(meta)[meta$`diabetes status:ch1` %in% c("ND", "T2D")]
res_all <- run_deseq(all_ids)

# B. Main sample set without the outlier GSM5009304
main_ids  <- rownames(meta)[meta$`diabetes status:ch1` %in% c("ND", "T2D") &
                              meta$`in_ins_filtered_data_subset:ch1` == "TRUE"]
res_noout <- run_deseq(setdiff(main_ids, "GSM5009304"))

sensitivity <- bind_rows(compare_to_main(res_all,   "All samples (n = 57)"),
                         compare_to_main(res_noout, "Without outlier (n = 38)"))
sensitivity
write.csv(sensitivity, "results/sensitivity_summary.csv", row.names = FALSE)