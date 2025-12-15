# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18")
# 
# BiocManager::install("BatchQC")

library(BatchQC)

pacman::p_load("dplyr", "readr", "stringr")

# Load raw count data
multicov_colnames.1 <- colnames(read.table("../results/Human_cells/D3_ATAC/D3_merged_peaks/celltype_condition.filteredNfixed.union.peakSet",
                                           sep = "\t", header = TRUE))[1:3]


multicov_colnames.2 <- read_tsv("../results/Human_cells/D3_ATAC/data/multicov_bam_order_master2.txt", col_names = "samples_id") %>% 
  mutate(samples_id = str_replace_all(samples_id, c("results/Human_cells/D3_ATAC/mapping/|_sorted_rmdup.bam" = ""))) %>% pull()

multicov_colnames <- c(multicov_colnames.1, multicov_colnames.2)

multicov <- read.table("../results/Human_cells/D3_ATAC/data/multicov_master_v2.tsv", sep = "\t",
                       col.names = c(multicov_colnames.1, multicov_colnames.2))



cts <- multicov %>% 
  mutate(peaks = paste(seqnames, start, end, sep = "-")) %>% 
  column_to_rownames("peaks") %>% 
  select(starts_with("WT_V"), starts_with("KD_V"),
         starts_with("WT_C"), starts_with("KD_C"),
         starts_with("W1_V"), starts_with("K1_V"),
         starts_with("W1_C"), starts_with("K1_C"),
         starts_with("WT_T"), starts_with("KD_T"),
         starts_with("WC_T"), starts_with("KC_T")) 


write.csv(cts, file = "../results/Human_cells/D3_ATAC/data/counts_ordered.csv")

exprdata <- as.matrix(cts)

metadf <- data.frame(Batch = c(rep(c(rep(c("A", "B"), c(2, 1))), 4), rep("C", 12), rep("B", 12) ), 
                     Group = rep(c(rep(c("WT", "KD"), c(3,3))), 6),
                     treatment = rep(c("veh", "ccm", "Veh1", "ccm1h", "t52", "ct52"), c(6, 6, 6, 6, 6, 6)),
                     row.names = colnames(exprdata))



table(metadf$Batch, interaction(metadf$Group, metadf$treatment, drop = TRUE))

# cols for AB cohort (veh/ccm only)
ab_cols <- rownames(metadf)[metadf$Batch %in% c("A","B") &
                              metadf$treatment %in% c("veh","ccm")]

cts_ab   <- as.matrix(cts[, ab_cols])
meta_ab  <- metadf[ab_cols, , drop = FALSE]

library(sva)
cb_counts_ab <- ComBat_seq(counts = cts_ab,
                           batch  = meta_ab$Batch,
                           group  = interaction(meta_ab$Group, meta_ab$treatment, drop = TRUE))


cts <- cbind(cb_counts_ab, exprdata[ ,13:36])
write.csv(cts, file = "../results/Human_cells/D3_ATAC/data/cts_batch_corrected.csv")

# Optionally save session info
sink("../logs/sessioninfo_human_09_batch_correction_ATAC_D3.txt")
sessionInfo()
sink()

