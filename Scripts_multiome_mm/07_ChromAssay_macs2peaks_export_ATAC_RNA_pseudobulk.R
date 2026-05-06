
# Load libraries
pacman::p_load("Seurat", "Signac", "harmony", "vroom", "ComplexHeatmap", "data.table",
               "SummarizedExperiment", "dplyr", "ggplot2", "rtracklayer", "RColorBrewer",
               "pheatmap", "BSgenome.Mmusculus.UCSC.mm10", "DoubletFinder", "EnsDb.Mmusculus.v79",
               "dittoSeq", "GenomicRanges", "colorspace", "readr", "stringr")

#### If starting with seurat object from GEO you can skip this part----

# # Load genome annotations
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# seqlevelsStyle(annotations) <- "UCSC"
# genome(annotations) <- "mm10"
# 
# # Load final Seurat object
# endo <- readRDS("../results/seurat_objects/merged/endo_annotated_final.rds")
# 
# # Load filtered union peak set
# peak_file <- "../results/ATAC/merged_peaks/celltype_condition.filteredNfixed.union.peakSet"
# peaks.df <- read.table(peak_file, sep = "\t", header = TRUE) %>%
#   dplyr::filter(spm >= 2)
# peaks <- makeGRangesFromDataFrame(peaks.df)
# 
# # Quantify fragments across peaks
# DefaultAssay(endo) <- "ATAC"
# macs2_counts <- FeatureMatrix(fragments = Fragments(endo),
#                               features = peaks,
#                               cells = colnames(endo),
#                               verbose = FALSE)
# 
# # Add quantified peaks as new assay
# endo[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts,
#                                         fragments = Fragments(endo),
#                                         annotation = annotations,
#                                         genome = "mm10",
#                                         verbose = FALSE)
# 
# DefaultAssay(endo) <- "peaks"
# 
# # Save Seurat object with links
# saveRDS(endo, "../results/seurat_objects/merged/endo_annotated_peaks_called.rds")
# 
# library(future)
# library(pryr)
# 
# tryCatch({
#   plan("multisession", workers = 2)
# }, error = function(e) {
#   message("⚠️ Parallel workers failed. Falling back to sequential.")
#   plan("sequential")
# })
# 
# options(future.globals.maxSize = 64000 * 1024^2)  # 64 GiB
# 
# # Compute peak–gene links
# # first compute the GC content for each peak
# DefaultAssay(endo) <- "peaks"
# endo <- RegionStats(endo, genome = BSgenome.Mmusculus.UCSC.mm10)
# 
# 
# start_time <- Sys.time()
# message("🔗 Running LinkPeaks...")
# endo <- LinkPeaks(object = endo,
#                   peak.assay = "peaks",
#                   expression.assay = "RNA",
#                   peak.slot = "counts",
#                   expression.slot = "data",
#                   pvalue_cutoff = 0.05)
# end_time <- Sys.time()
# duration <- end_time - start_time
# message(glue::glue("✅ LinkPeaks completed in {round(duration, 2)}"))
# message(glue::glue("🔍 Memory used: {round(pryr::mem_used() / 1024^3, 2)} GB"))
# 
# # Export peak–gene links
# links <- endo@assays$peaks@links
# links_df <- as.data.frame(links)
# colnames(links_df)[colnames(links_df) == "peak"] <- "peaks"
# write.csv(links_df, "../results/ATAC/data/endo_linkpeaks.csv", row.names = FALSE)
# 
# # Save Seurat object with links
# saveRDS(endo, "../results/seurat_objects/merged/endo_annot_MACS2peaks_wLinks.rds")

### Starting with RDS object from GEO ----
endo <- readRDS("../results/seurat_objects/merged/endo_annot_MACS2peaks_wLinks.rds")


# Aggregate ATAC (pseudobulk) by cell type and sample----
endo$samples <- paste(endo$cell_type, endo$orig.ident, sep = "_")
cts_atac <- AggregateExpression(endo,
                                group.by = "samples",
                                assays = "peaks",
                                return.seurat = FALSE)$peaks

cts_atac <- as.data.frame(cts_atac)
colnames(cts_atac) <- gsub(" ", "_", colnames(cts_atac))
cts_atac <- cts_atac[c(1:12, 19:24, 13:18, 25:30)]


# Export ATAC counts
dir.create("../results/ATAC/data", recursive = TRUE, showWarnings = FALSE)
write.csv(cts_atac, "../results/ATAC/data/counts_atac_all.csv", row.names = TRUE)


# Aggregate RNA pseudobulk----
DefaultAssay(endo) <- "RNA"
cts_rna <- AggregateExpression(endo,
                               group.by = "samples",
                               assays = "RNA",
                               return.seurat = FALSE)$RNA

cts_rna <- as.data.frame(cts_rna)
colnames(cts_rna) <- gsub(" ", "_", colnames(cts_rna))
cts_rna <- cts_rna[c(1:12, 19:24, 13:18, 25:30)]

# Export RNA counts
dir.create("../results/RNA/data", recursive = TRUE, showWarnings = FALSE)
write.csv(cts_rna, "../results/RNA/data/counts_rna_all.csv", row.names = TRUE)


annot_peaks <- read_tsv("../results/ATAC/data/annotated_union_clean_sorted.peakset.txt") %>% 
  rename(peaks = starts_with("PeakID")) %>% 
  dplyr::select(peaks, Annotation) %>% 
  mutate(Annotation_final = case_when(str_detect(Annotation, regex("promoter", ignore_case = TRUE)) ~ "Promoter-TSS",
                                      str_detect(Annotation, regex("5' UTR", ignore_case = TRUE)) ~ "5' UTR",
                                      str_detect(Annotation, regex("3' UTR", ignore_case = TRUE)) ~ "3' UTR",
                                      str_detect(Annotation, regex("exon", ignore_case = TRUE)) ~ "Exon",
                                      str_detect(Annotation, regex("Intron", ignore_case = TRUE)) ~ "Intron",
                                      str_detect(Annotation, regex("intergenic", ignore_case = TRUE)) ~ "Distal\nIntergenic",
                                      str_detect(Annotation, regex("TTS", ignore_case = TRUE)) ~ "TTS",
                                      TRUE ~ Annotation))


annotation_counts <- annot_peaks %>% 
  group_by(Annotation_final) %>% 
  summarise(Count = n()) %>% 
  mutate(Percentage = round((Count / sum(Count)) * 100, 1))

annotation_counts$Annotation_final <- factor(annotation_counts$Annotation_final, levels = c("Promoter-TSS", "5' UTR", "3' UTR", "Exon", 
                                                                                            "Intron", "Distal\nIntergenic", "TTS"))
annotation_counts <- annotation_counts %>% 
  arrange(Annotation_final)

# Generate the color palette
color_mapping <- qualitative_hcl(7, palette = "Set 3")

annotation_counts <- annotation_counts %>%
  mutate(label = paste(Annotation_final, " (", round(Percentage, 1), "%)", sep = ""))
annotation_counts$label <- factor(annotation_counts$label, levels = unique(annotation_counts$label))


# Adjust the color mapping to the new labels
color_mapping <- setNames(color_mapping, annotation_counts$label)

annot_peaks_flox_enhancers <- ggplot(annotation_counts, aes(x = "", y = Percentage, fill = label)) +
  geom_bar(width = 1, stat = "identity", color = "white" ) +
  coord_polar("y") + 
  scale_fill_manual(values = color_mapping) +
  theme_void() + 
  labs(title = "Genomic annotation of non-redundant BEC cCRE atlas",
       caption = paste0(sum(annotation_counts$Count), " cis-regulatory elements (cCRE)")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        legend.title = element_blank(), plot.caption.position = "plot", plot.caption = element_text(size = 9))

ggsave(filename = "../results/plots/graphs_ATAC/non_redundant_cCRE_annotation.svg", plot = annot_peaks_flox_enhancers, 
       device = "svg", width = 4.8, height = 2.5, units = "in", bg = "white")


# Save session info
dir.create("../results/logs", showWarnings = FALSE)
sink("../results/logs/sessioninfo_07_ChromAssay_macs2peaks_export_ATAC_RNA_pseudobulk.txt")
sessionInfo()
sink()
