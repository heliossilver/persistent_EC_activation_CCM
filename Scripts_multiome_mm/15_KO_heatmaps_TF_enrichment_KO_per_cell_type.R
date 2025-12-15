# Load required libraries
pacman::p_load("tibble", "RColorBrewer", "pheatmap", "readr",
               "stringr", "purrr", "viridis", "vroom", "dendsort",
               "ComplexHeatmap", "circlize", "glue", "patchwork",
               "dplyr", "ggplot2", "textshape", "paletteer", "tidyr")

#### Step 1: Clean HOMER output ----

motif_dir <- "../results/RNA_ATAC_integration/Homer/enhancers_KO_enrichment/"
folders <- list.dirs(motif_dir, full.names = FALSE, recursive = FALSE)
folders <- folders[c(1, 3, 2, 4, 5)]  # Reorder manually

for (folder in folders) {
  
  motif_file <- glue("{motif_dir}/{folder}/knownResults.txt")
  new_name <- gsub("_homer", "", folder)
  
  motif_df <- vroom::vroom(motif_file) %>%
    select(`Motif Name`, `P-value`, `% of Target Sequences with Motif`, `% of Background Sequences with Motif`) %>%
    rename(TF_motif = `Motif Name`,
           P_value = `P-value`) %>%
    mutate(log10_pvalue = -log10(as.numeric(P_value)),
           TF_motif = str_extract(TF_motif, "^[^(]+"))
  
  write.csv(motif_df, file = glue("{motif_dir}/{folder}/{new_name}.csv"), row.names = FALSE)
}

#### Step 2: Load and filter cleaned motif data ----

samples <- gsub("_homer", "", folders)
threshold <- -log10(0.05)

motif_list <- map(samples, function(sample) {
  read_csv(glue("{motif_dir}/{sample}_homer/{sample}.csv")) %>%
    mutate(cell_type = sample) %>%
    filter(log10_pvalue > threshold) %>%
    select(TF_motif, log10_pvalue, cell_type) %>%
    distinct(TF_motif, .keep_all = TRUE)
})

names(motif_list) <- samples

motif_df <- bind_rows(motif_list)
motif_df$cell_type <- factor(motif_df$cell_type, levels = samples)

write.csv(motif_df, file = "../results/RNA_ATAC_integration/Homer/results_tables/enhancer_ko_cCRE_knownTF_percell.csv", row.names = FALSE)

#### Step 3: Full heatmap ----

motif_matrix <- motif_df %>%
  pivot_wider(id_cols = TF_motif, names_from = cell_type, values_from = log10_pvalue, values_fill = 0) %>%
  column_to_rownames("TF_motif")

order_motifs <- motif_df %>%
  group_by(TF_motif) %>%
  arrange(desc(log10_pvalue), .by_group = TRUE) %>%
  distinct(TF_motif, .keep_all = TRUE) %>%
  ungroup() %>%
  arrange(cell_type, desc(log10_pvalue)) %>%
  pull("TF_motif")

motif_matrix <- motif_matrix[order_motifs, ]

annot_df <- data.frame(`Endothelial cell` = c("Artery KO", "Capillary artery KO", "Capillary KO", "Capillary vein KO", "Large Vein KO"),
  row.names = colnames(motif_matrix))

annot_df$Endothelial.cell<- factor(annot_df$Endothelial.cell, levels = unique(annot_df$Endothelial.cell))

EC_colors <- c("#247567", "#E2705B", "#6C7A99", "#B4698B", "#D8A21F")
names(EC_colors) <- levels(annot_df$Endothelial.cell)
annot_colors <- list("Endothelial cell" = EC_colors)

top_annot <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = "#FFF850", col = "#FFDA03")),
                               `Endothelial cell` = annot_df$Endothelial.cell,
                               col = list(`Endothelial cell` = EC_colors),
                               annotation_height = unit(c(2.5, 2.5), "mm"),
                               show_legend = c(FALSE, FALSE),
                               show_annotation_name = FALSE)

row_dend_full <- dendsort::dendsort(hclust(dist(motif_matrix)),  isReverse = TRUE)
hmap_full <- ComplexHeatmap::pheatmap(mat = motif_matrix,
                                      name = "-log10\nP-value",
                                      cluster_rows = row_dend_full,
                                      cluster_cols = FALSE,
                                      show_row_dend = FALSE,
                                      show_colnames = FALSE,
                                      show_rownames = TRUE,
                                      top_annotation = top_annot,
                                      heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                                  labels_gp = gpar(fontsize = 7),
                                                                  legend_height = unit(1.8, "cm")),
                                      color = c("white", rocket(150, begin = 1, end = 0)[1:99]),
                                      border = "black",
                                      annotation_colors = annot_colors,
                                      border_color = NA,
                                      main = "Transcription factor\nmotif enrichment")

png(filename = "../results/RNA_ATAC_integration/graphs/hmap_ko_known_names.png", width = 4, height = 24, units = "in", res = 300)
ComplexHeatmap::draw(hmap_full, heatmap_legend_side = "left", show_annotation_legend = FALSE)
dev.off()

#### Step 4: Top 15 motifs per cell type ----

motif_top15 <- motif_df %>%
  group_by(cell_type) %>%
  slice_max(order_by = log10_pvalue, n = 15) %>%
  ungroup() %>%
  group_by(TF_motif) %>%
  arrange(desc(log10_pvalue), .by_group = TRUE) %>%
  distinct(TF_motif, .keep_all = TRUE) %>%
  ungroup() %>%
  arrange(cell_type, desc(log10_pvalue)) %>%
  pull("TF_motif")

# Optionally append known TFs manually
motif_top15 <- unique(c(motif_top15, "NFkB-p65", "NFkB-p65-Rel", "NFkB2-p52"))

motif_matrix_top15 <- motif_df %>%
  pivot_wider(id_cols = TF_motif, names_from = cell_type, values_from = log10_pvalue, values_fill = 0) %>%
  filter(TF_motif %in% motif_top15) %>%
  column_to_rownames("TF_motif")

row_dend_15 <- dendsort::dendsort(hclust(dist(motif_matrix_top15)))

hmap_top15 <- ComplexHeatmap::pheatmap(mat = motif_matrix_top15,
                                       name = "-log10\nP-value",
                                       cluster_rows = row_dend_15,
                                       cluster_cols = FALSE,
                                       show_colnames = FALSE,
                                       show_rownames = TRUE,
                                       show_row_dend = FALSE,
                                       top_annotation = top_annot,
                                       heatmap_legend_param = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
                                                                   labels_gp = gpar(fontsize = 5),
                                                                   legend_height = unit(1.5, "cm")),
                                       breaks = seq(0, max(motif_matrix_top15), length.out = 101),
                                       color = c("white", rocket(150, begin = 1, end = 0)[1:99]),
                                       border = "black",
                                       border_color = "#FFFFFF",
                                       fontsize = 7,
                                       fontsize_row = 6,
                                       annotation_colors = annot_colors,
                                       main = "Transcription factor\n motif enrichment")

png(filename = "../results/RNA_ATAC_integration/graphs/hmap_ko_top15_known.png", width = 2.5, height = 3, units = "in", res = 300)
ComplexHeatmap::draw(hmap_top15, merge_legends = TRUE, heatmap_legend_side = "left", show_annotation_legend = FALSE)
dev.off()

# === Log R session info for reproducibility ===
sink("../logs/sessioninfo_15_KO_heatmaps_TF_enrichment_KO_per_cell_type.txt")
sessionInfo()
sink()

