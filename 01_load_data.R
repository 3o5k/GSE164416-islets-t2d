# 01_load_data.R
# Load NCBI-generated raw counts and GEO sample metadata for GSE164416 and
# define the analysis set: ND vs T2D, samples passing the INS purity filter.
# Input : data/GSE164416_raw_counts_GRCh38.p13_NCBI.tsv.gz (downloaded if missing)
# Output: results/counts_sub.rds, results/coldata.rds, results/meta.rds

library(tidyverse)
library(GEOquery)

dir.create("data",    showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# 0. Download the NCBI-generated count matrix and gene annotation if missing
counts_file <- "data/GSE164416_raw_counts_GRCh38.p13_NCBI.tsv.gz"
annot_file  <- "data/Human.GRCh38.p13.annot.tsv.gz"

if (!file.exists(counts_file)) {
  download.file(paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts",
                       "&acc=GSE164416&format=file&file=", basename(counts_file)),
                counts_file, mode = "wb")
}
if (!file.exists(annot_file)) {
  download.file(paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts",
                       "&format=file&file=", basename(annot_file)),
                annot_file, mode = "wb")
}

# 1. Count matrix: rows = genes (Entrez IDs), columns = samples (GSM IDs)
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
dim(counts)
counts[1:5, 1:5]

# 2. Sample metadata from GEO
gse  <- getGEO("GSE164416", GSEMatrix = TRUE)
meta <- pData(gse[[1]])
dim(meta)

# 3. Inspect the key metadata columns
table(meta$`diabetes status:ch1`)
table(meta$`in_ins_filtered_data_subset:ch1`)
table(meta$`library_prep_date:ch1`)

# 4. Check that sample order matches between counts and metadata
all(colnames(counts) == meta$geo_accession)

# 5. Cross-tabulate diagnosis against batch and the INS-filter flag
table(meta$`diabetes status:ch1`, meta$`library_prep_date:ch1`)
table(meta$`diabetes status:ch1`, meta$`in_ins_filtered_data_subset:ch1`)

# 6. Does the INS-filter flag reflect insulin expression? (INS Entrez ID = 3630)
ins_cpm <- as.numeric(counts["3630", ]) / colSums(counts) * 1e6
tapply(ins_cpm, meta$`in_ins_filtered_data_subset:ch1`, median)

# 7. Analysis subset: ND vs T2D, INS-filtered (purity-passing) samples
keep <- meta$`diabetes status:ch1` %in% c("ND", "T2D") &
  meta$`in_ins_filtered_data_subset:ch1` == "TRUE"
sum(keep)

counts_sub <- counts[, keep]

coldata <- data.frame(
  condition = factor(meta$`diabetes status:ch1`[keep], levels = c("ND", "T2D")),
  batch     = factor(meta$`library_prep_date:ch1`[keep]),
  row.names = meta$geo_accession[keep]
)

dim(counts_sub)
table(coldata$condition, coldata$batch)
all(colnames(counts_sub) == rownames(coldata))

# 8. Save the analysis objects for the next scripts
saveRDS(counts_sub, "results/counts_sub.rds")
saveRDS(coldata,    "results/coldata.rds")
saveRDS(meta,       "results/meta.rds")