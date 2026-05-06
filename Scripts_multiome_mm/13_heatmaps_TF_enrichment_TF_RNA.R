
# Load required packages
pacman::p_load("tibble", "RColorBrewer", "stringr", "purrr", "viridis", "ggtext",
               "ComplexHeatmap", "circlize", "glue", "patchwork", "readr",
               "dplyr", "ggplot2", "textshape", "paletteer", "tidyr")

#### 1. Clean motif enrichment output files from HOMER ----
# Define input and output paths
motif_dir <- "../results/RNA_ATAC_integration/Homer/enhancer_flox_enrichment"
output_csv_dir <- motif_dir  # Cleaned CSVs will be stored alongside raw files

results_tables_homer_dir <- "../results/RNA_ATAC_integration/Homer/results_tables"
dir.create(results_tables_homer_dir, recursive = TRUE, showWarnings = FALSE)


# List cell-type specific motif result folders
folders <- list.dirs(motif_dir, full.names = FALSE, recursive = FALSE)
folders <- folders[c(1, 3, 2, 4, 5)]  # Reorder manually if needed

# Loop through each file and clean the knownResults.txt
for (folder in folders) {
  
  motif_file <- glue("{motif_dir}/{folder}/knownResults.txt")
  new_name <- gsub("_markers_homer", "", folder)
  
  motif_df <- vroom::vroom(motif_file) %>% 
    select(`Motif Name`, `P-value`, `% of Target Sequences with Motif`) %>% 
    rename(TF_motif = `Motif Name`,
      P_value = `P-value`,
      Percent_Targets = `% of Target Sequences with Motif`) %>% 
    mutate(log10_pvalue = -log10(as.numeric(P_value)),
      Percent_Targets = as.numeric(gsub("%", "", Percent_Targets)),
      TF_motif = str_extract(TF_motif, "^[^(]+"))  # Remove trailing details
  
  # Write cleaned CSV
  write.csv(motif_df, file = glue("{motif_dir}/{folder}/{new_name}.csv"), row.names = FALSE)
}


#### 2. Load cleaned motif data and prepare for visualization ----

samples <- gsub("_markers_homer", "", folders)

# Precalculate -log10(0.05) threshold
log10_thresh <- -log10(0.05)

# Read and filter cleaned motif tables
motif_list <- map(samples, function(sample) {
  read_csv(glue("{motif_dir}/{sample}_markers_homer/{sample}.csv")) %>%
    filter(log10_pvalue > log10_thresh) %>%
    mutate(cell_type = sample) %>%
    distinct(TF_motif, .keep_all = TRUE) %>%
    select(TF_motif, log10_pvalue, cell_type)
})

names(motif_list) <- samples

# Combine into one data frame
motif_long <- bind_rows(motif_list)
motif_long$cell_type <- factor(motif_long$cell_type, levels = samples)

# Determine ordering of TF motifs for heatmap
motif_order <- motif_long %>%
  group_by(TF_motif) %>%
  arrange(desc(log10_pvalue), .by_group = TRUE) %>%
  distinct(TF_motif, .keep_all = TRUE) %>%
  ungroup() %>%
  arrange(cell_type, desc(log10_pvalue)) %>%
  pull(TF_motif)

# Convert to wide format for heatmap
motif_wide <- motif_long %>%
  pivot_wider(names_from = cell_type, values_from = log10_pvalue, values_fill = 0) %>%
  column_to_rownames("TF_motif")

# Calculate row gaps for heatmap grouping
row_gaps <- motif_long %>%
  group_by(TF_motif) %>%
  arrange(desc(log10_pvalue), .by_group = TRUE) %>%
  distinct(TF_motif, .keep_all = TRUE) %>%
  ungroup() %>%
  arrange(cell_type, desc(log10_pvalue)) %>%
  group_by(cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(cumsum_n = cumsum(n)) %>%
  pull(cumsum_n)


# Save long-format table to data folder
write.csv(motif_long, file = glue("{results_tables_homer_dir}/enhancer_flox_known_TF_percell.csv"), row.names = FALSE)


#### RNA TF lvls ----
#getting orthologs from jax


jax_ortholog_path <- "utils/jax_ms_hu_orthologs.tsv"
if (!file.exists(jax_ortholog_path)) {
  download.file("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
                destfile = jax_ortholog_path)
}
mouse_human_genes <- read.table(jax_ortholog_path, sep = "\t", header = TRUE)

convert_human_to_mouse <- function(gene_list){
  human_filtered <- mouse_human_genes %>% filter(Symbol %in% gene_list, Common.Organism.Name == "human" )
  mouse_filtered <- mouse_human_genes %>% filter(Common.Organism.Name == "mouse, laboratory" )
  merged <- left_join(human_filtered, mouse_filtered, by = "DB.Class.Key") %>%
    dplyr::select(Symbol.x, Symbol.y) %>%
    rename(gene_hu = Symbol.x, gene_ms = Symbol.y)
  return(merged)
}

# Split motifs that include dimerization (e.g., "Fos-Jun")
motifs_split <- motif_long %>%
  separate_rows(TF_motif, sep = "[-:]")  # Handles "-" and ":" split

motifs <- pull(motifs_split, "TF_motif")

# Categorize motif names
human_motifs  <- motifs[grepl("^[A-Z0-9]+$", motifs)]
mouse_motifs  <- motifs[grepl("^[A-Z][a-z0-9]+$", motifs)]
other_motifs  <- motifs[!motifs %in% c(human_motifs, mouse_motifs)]

mouse_genes <- mouse_human_genes %>%
  filter(Common.Organism.Name == "mouse, laboratory") %>% pull(Symbol)

human_motifs_df <- data.frame(gene_hu = human_motifs)
human_to_mouse_motifs <- convert_human_to_mouse(human_motifs)

human_motifs_df <- human_motifs_df %>%
  left_join(human_to_mouse_motifs, by = "gene_hu")

# Collect matched TFs
motifs_mouse_tolook <- mouse_motifs[mouse_motifs %in% mouse_genes]
motifs_converted <- human_motifs_df %>% drop_na() %>% pull(gene_ms)
motifs_mouse_tolook <- unique(c(motifs_mouse_tolook, motifs_converted))

# Capture unmatched motifs for manual matching
motifs_not_matched_ms <- setdiff(mouse_motifs, mouse_genes)
motif_not_matched_hu <- human_motifs_df %>%
  filter(is.na(gene_ms)) %>% pull(gene_hu)

# Add likely matches for ambiguous motifs
manual_matches <- c("Fosl1", "Fosl2", "Nfe2l2", "Nfe2", "Rbpj", "Tfap4", "Nr4a1", "Esrra",
                    "Erg", "Fev", "Ear2", "Nr2c2", "Ctcfl", "Cebpb", "Tal1", "Pgr", "Mtf1")


# Add gene families from RNA data
files_rna <- list.files("../results/RNA/edgeR/deg_tables/", pattern = "_vs_rest.csv")
files_rna <- files_rna[c(1, 3, 2, 4, 5)]
files_rna <- gsub(".csv", "", files_rna)

datalist <- map(files_rna, function(f) {
  vroom::vroom(file = glue("../results/RNA/edgeR/deg_tables/{f}.csv")) %>%
    filter(get(glue("FDR_{f}")) <= 0.05,
           abs(get(glue("logFC_{f}"))) > 0) %>%
    arrange(get(glue("FDR_{f}"))) %>%
    select(gene, starts_with("logFC"))
})
names(datalist) <- files_rna

# Merge all RNA logFCs
merged_df <- reduce(datalist, full_join, by = "gene") %>%
  replace(is.na(.), 0)

# Additional TF gene families (from merged_df)
family_genes <- c(grep("^Oct", merged_df$gene, value = TRUE),
                  grep("^Fox", merged_df$gene, value = TRUE),
                  grep("^Zic", merged_df$gene, value = TRUE),
                  grep("^Klf", merged_df$gene, value = TRUE),
                  grep("^Ebf", merged_df$gene, value = TRUE),
                  grep("^Ets", merged_df$gene, value = TRUE),
                  grep("^Etv", merged_df$gene, value = TRUE),
                  grep("^Nfat", merged_df$gene, value = TRUE),
                  grep("^Gata", merged_df$gene, value = TRUE),
                  grep("^Tead", merged_df$gene, value = TRUE),
                  grep("^Rfx", merged_df$gene, value = TRUE),
                  grep("^Irf", merged_df$gene, value = TRUE),
                  grep("^Tcf", merged_df$gene, value = TRUE),
                  grep("^Stat", merged_df$gene, value = TRUE),
                  grep("^Rxr", merged_df$gene, value = TRUE),
                  grep("^Runx", merged_df$gene, value = TRUE),
                  grep("^Nfkb", merged_df$gene, value = TRUE))


ap1_genes <- c("Jun", "Junb", "Jund", "Fos", "Fosb", "Atf2", "Atf3", "Batf", "Dnajc12",
               "dp2", "Maf", "Mafb", "Mafa", "Mafg" ,"Maff", "Mafk", "Nrl")

# Final list of genes to look at
motifs_mouse_tolook <- unique(c(motifs_mouse_tolook, manual_matches, family_genes, ap1_genes))

# Filter for these TFs
merged_df_filt <- merged_df %>%
  filter(gene %in% motifs_mouse_tolook) %>%
  column_to_rownames("gene") %>%
  filter(if_any(everything(), ~ . > 0))

# Long format for ordering
merged_df_filt_long <- merged_df_filt %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell_type", values_to = "logFC") %>%
  drop_na() %>%
  group_by(gene) %>%
  arrange(desc(logFC), .by_group = TRUE) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(logFC))

merged_df_filt_long$cell_type <- factor(merged_df_filt_long$cell_type, levels = colnames(merged_df_filt))

# Order for heatmap
gene_order <- merged_df_filt_long %>%
  arrange(cell_type) %>%
  pull(gene)

# Save final expression table of predicted TFs
write.csv(merged_df_filt[gene_order, ], file = glue("{results_tables_homer_dir}/RNA_FC_knownTF.csv"), row.names = TRUE)


###HEATMAPS
#### 1. Get labels for TFs heatmaps----

merged_df_filt_long_clean <- merged_df_filt_long %>% 
  arrange(desc(logFC)) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  mutate(gene_lower = tolower(gene)) %>% select(-cell_type)

human_to_mouse_motifs <- bind_rows(human_to_mouse_motifs, data.frame(gene_hu = c("NFkB-p65", "NFkB-p65-Rel"),
                                                                     gene_ms = c("Nfkb1", "Rel")))

motif_long_clean <- motif_long %>% 
  arrange(desc(log10_pvalue)) %>% 
  distinct(TF_motif, .keep_all = TRUE) %>%
  left_join(human_to_mouse_motifs, by = c("TF_motif" = "gene_hu")) %>%
  mutate(TF_motif_to_join = coalesce(gene_ms, TF_motif)) %>%
  select(-gene_ms) %>% 
  mutate(TF_motif_to_join = tolower(TF_motif_to_join)) 


# Merge by gene/motif and cell_type
# TF_motif from motif enrichment should match gene from RNA
labels_TF_RNA <- inner_join(motif_long_clean,
                        merged_df_filt_long_clean,
                        by = c("TF_motif_to_join" = "gene_lower")) %>% 
  select(TF_motif, gene, cell_type, log10_pvalue) %>%
  rename(TF_RNA = gene) %>%
  distinct() %>% 
  arrange(cell_type, desc(log10_pvalue)) %>% 
  group_by(cell_type) 


labels_to_keep <- c(labels_TF_RNA %>%  
  slice_head(n = 7) %>% pull(TF_RNA), "Gata6", "Stat4", "Tcf19", "Rel") 

labels_to_remove <- c("Klf6", "Klf3", "Klf10", "Zbtb18")

labels_df <- labels_TF_RNA %>% 
  filter(TF_RNA %in% labels_to_keep,
         !TF_RNA %in% labels_to_remove)


#### 2. Motif enrichment heatmap ----

# Annotation: Endothelial cell types
EC_list <- c("Artery", "Capillary artery", "Capillary", "Capillary vein", "Large Vein")
annot_tf <- data.frame(`Endothelial cell` = EC_list, row.names = colnames(motif_wide))

annot_tf$Endothelial.cell <- factor(annot_tf$Endothelial.cell, levels = EC_list)

# Colors
EC_colors <- c("#247567", "#E2705B", "#6C7A99", "#B4698B", "#D8A21F")
names(EC_colors) <- EC_list
my_colour <- list("Endothelial cell" = EC_colors)

# Order matrix by TFs
ec_motif_full_wide <- motif_wide[motif_order, ]
rownames_ec_motif <- rownames(ec_motif_full_wide)

# Label positions
tf_positions <- match(labels_df$TF_motif, rownames_ec_motif)
tf_positions_anno <- tf_positions[!is.na(tf_positions)]

# Right-side annotation
ha <- rowAnnotation(motif_labels = anno_mark(at = tf_positions_anno,
                                             labels = labels_df$TF_motif[!is.na(tf_positions)],
                                             labels_gp = gpar(fontsize = 7.5, col = "black"),
                                             padding = unit(1.2, "mm")))

# Top annotation (cell type blocks)
top_annot_tf <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = "#FFF850", col = "#FFDA03")),
                                  `Endothelial cell` = annot_tf$Endothelial.cell,
                                  col = list(`Endothelial cell` = EC_colors),
                                  annotation_height = unit(c(2.5, 2.5), "mm"),
                                  show_legend = c(FALSE, FALSE),
                                  show_annotation_name = FALSE)

# Heatmap
hmap_flox_markers <- ComplexHeatmap::pheatmap(mat = ec_motif_full_wide,
                                              name = "-log10\nP-value",
                                              cluster_rows = FALSE,
                                              cluster_cols = FALSE,
                                              show_colnames = FALSE,
                                              show_rownames = FALSE,
                                              heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                                          labels_gp = gpar(fontsize = 7),
                                                                          legend_height = unit(1.8, "cm")),
                                              top_annotation = top_annot_tf,
                                              right_annotation = ha,
                                              breaks = seq(-1, max(ec_motif_full_wide), length.out = 101),
                                              #color = colorRampPalette(brewer.pal(9, "Reds"))(100),
                                              color = c("white", rocket(150, begin = 1, end = 0)[1:99]),
                                              main = "Transcription factor\nmotif enrichment",
                                              border_color = NA,
                                              border = "black",
                                              annotation_names_col = FALSE,
                                              fontsize = 8,
                                              fontsize_row = 6,
                                              gaps_row = row_gaps,
                                              annotation_colors = my_colour)

# Save heatmap to file
png(filename = "../results/RNA_ATAC_integration/graphs/hmap_flox_known.png", width = 2.5, height = 4, units = "in", res = 300)
ComplexHeatmap::draw(hmap_flox_markers, heatmap_legend_side = "left", show_annotation_legend = FALSE)
dev.off()

#### 3. RNA logFC heatmap ----

annot_gTF <- data.frame(`Endothelial cell` = EC_list, row.names = colnames(merged_df_filt))
annot_gTF$Endothelial.cell <- factor(annot_gTF$Endothelial.cell, levels = EC_list)

# Match rownames for expression matrix
merged_df_filt_wide <- merged_df_filt[gene_order, ]
rownames_ec_merged_df_filt <- rownames(merged_df_filt_wide)

# Label positions
RNA_labels <- c(labels_df$TF_RNA, "Tcf19")
gTF_positions <- match(RNA_labels, rownames_ec_merged_df_filt)
gTF_positions_anno <- gTF_positions[!is.na(gTF_positions)]

# Right annotation for RNA heatmap
ha_gTF <- rowAnnotation(motif_labels = anno_mark(at = gTF_positions_anno,
                                                 labels = RNA_labels[!is.na(gTF_positions)],
                                                 labels_gp = gpar(fontsize = 7.5, col = "black"),
                                                 padding = unit(1.2, "mm")))

# Top annotation
top_annot_gTF <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = "#BEE880", col = "#027000")),
                                   `Endothelial cell` = annot_gTF$Endothelial.cell,
                                   col = list(`Endothelial cell` = EC_colors),
                                   annotation_height = unit(c(2.5, 2.5), "mm"),
                                   show_legend = c(FALSE, FALSE),
                                   show_annotation_name = FALSE)

# Row gap locations based on long format expression
gaps_row_gTF_15 <- merged_df_filt_long %>%
  group_by(cell_type) %>%
  summarise(n = n()) %>%
  mutate(cum = cumsum(n)) %>%
  pull(cum) %>%
  head(4)

# Heatmap for RNA expression of TFs
hmap_flox_TF_RNA <- ComplexHeatmap::pheatmap(mat = merged_df_filt_wide,
                                             name = "logFC",
                                             cluster_rows = FALSE,
                                             cluster_cols = FALSE,
                                             show_colnames = FALSE,
                                             show_rownames = FALSE,
                                             top_annotation = top_annot_gTF,
                                             right_annotation = ha_gTF,
                                             heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                                         labels_gp = gpar(fontsize = 7),
                                                                         legend_height = unit(1.8, "cm")),
                                             color = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)),
                                             breaks = seq(-3.5, 3.5, length.out = 101),
                                             main = "Transcription factors\ngene expression",
                                             border_color = "#FFFFFF",
                                             border = "black",
                                             annotation_colors = my_colour,
                                             fontsize = 8,
                                             fontsize_row = 6,
                                             gaps_row = gaps_row_gTF_15)

# Save RNA heatmap to file
png(filename = "../results/RNA_ATAC_integration/graphs/hmap_flox_TF_RNA_known.png", width = 2.5, height = 4, units = "in", res = 300)
ComplexHeatmap::draw(hmap_flox_TF_RNA, merge_legends = TRUE, heatmap_legend_side = "left", show_annotation_legend = FALSE)
dev.off()



#combined image ----
# Read images back and combine them horizontally
img_gene <- magick::image_read("../results/RNA_ATAC_integration/graphs/flox_cCRE_gene_2hmaps.png")
img_tf <- magick::image_read("../results/RNA_ATAC_integration/graphs/hmap_flox_known.png")
img_gTF <- magick::image_read("../results/RNA_ATAC_integration/graphs/hmap_flox_TF_RNA_known.png")

combined <- magick::image_append(c(img_gene, img_tf, img_gTF), stack = FALSE)
magick::image_write(combined, path = "../results/RNA_ATAC_integration/graphs/flox_heatmaps_combined_known.svg", format = "svg")


# === Log R session info for reproducibility ===
sink("../results/logs/sessioninfo_13_heatmaps_TF_enrichment_TF_RNA.txt")
sessionInfo()
sink()
