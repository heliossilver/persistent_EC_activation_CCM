
pacman::p_load("tibble", "RColorBrewer", "pheatmap",
               "stringr", "glue", "purrr", "dplyr",
               "ggplot2", "textshape", "patchwork",
               "paletteer", "tidyr","circlize", 
               "viridis", "ComplexHeatmap", "grid", "gridExtra")

# ==== Load and scale log2CPM values ====

log2cpm_atac <- read.csv(file = "../results/ATAC/data/log2cpm.csv") %>% column_to_rownames("peaks") %>% 
  dplyr::select(contains("Flox"))
log2cpm.scaled_atac <- as.data.frame(t(apply(log2cpm_atac, 1, scale)))

colnames(log2cpm.scaled_atac) <- colnames(log2cpm_atac)
log2cpm.scaled_atac <- rownames_to_column(.data = log2cpm.scaled_atac, "peaks")


log2cpm_rna <- read.csv(file = "../results/RNA/data/log2cpm.csv") %>% column_to_rownames("gene") %>% 
  dplyr::select(contains("Flox"))
log2cpm.scaled_rna <- as.data.frame(t(apply(log2cpm_rna, 1, scale)))

colnames(log2cpm.scaled_rna) <- colnames(log2cpm_rna)
log2cpm.scaled_rna <- rownames_to_column(.data = log2cpm.scaled_rna, "gene")
genes.matrix_flox <- log2cpm.scaled_rna$gene

# Peaks used in DARs
files <- list.files(path = "../results/ATAC/edgeR/dar_diffacc_tables", pattern = "_vs_rest_diffacc.csv")
samples <- gsub("_vs_rest_diffacc.csv", "", files[c(1, 3, 2, 4, 5)])
peaks.v <- c()
for (sample in samples) {
  data <- read.csv(file = sprintf("../results/ATAC/edgeR/dar_diffacc_tables/%s_vs_rest_diffacc.csv", sample)) %>%
    pull(peaks)
  peaks.v <- unique(c(peaks.v, data))
}


# Load pathway gene sets
relevant_genes <- read.csv("utils/relevante_genes_mouse.csv", fill = TRUE, na.strings = "") %>%
  select(-starts_with("endoMT"), -contains("old"))

# Filter genes function
filter_genes_f <- function(genes_df, pathway, genes_in_matrix) {
  selected_genes <- genes_df %>%
    select({{pathway}}) %>% drop_na() %>%
    filter({{pathway}} %in% genes_in_matrix) %>%
    pull({{pathway}})
  return(selected_genes)
}

# Match genes to matrix and save
pathway_genes <- map(names(relevant_genes), ~ filter_genes_f(relevant_genes, !!sym(.x), genes.matrix_flox))
names(pathway_genes) <- colnames(relevant_genes)
saveRDS(pathway_genes, file = "../results/RNA/rds/pathway_genes_inflox.rds")

# Load and filter DEGs for each flox sample
deg_files <- list.files("../results/RNA/edgeR/deg_tables", pattern = "_vs_rest.csv", full.names = TRUE)
deg_files <- deg_files[c(1, 3, 2, 4, 5)]  # reordered manually

samples <- gsub(".csv", "", basename(deg_files))
datalist <- setNames(lapply(samples, function(sample) {
  read.csv(glue("../results/RNA/edgeR/deg_tables/{sample}.csv")) %>%
    filter(get(glue("FDR_{sample}")) <= 0.05,
           get(glue("logFC_{sample}")) >= 0.5) %>%
    select(gene, starts_with("logFC")) %>%
    rename(logFC = starts_with("logFC")) %>%
    arrange(desc(logFC)) %>%
    mutate(cell_type = gsub("_vs_rest", "", sample))
}), samples)

# Ordered list and wide logFC tables per pathway
custom_order <- c("Artery", "CapArt", "Cap", "CapVein", "LargeVein")

get_ordered_genes <- function(pathway_genes, datalist, custom_order) {
  bind_rows(lapply(names(datalist), function(ct) {
    datalist[[ct]] %>%
      filter(gene %in% pathway_genes) %>%
      mutate(cell_type = factor(cell_type, levels = custom_order))
  })) %>%
    arrange(desc(logFC)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    arrange(cell_type) %>%
    select(gene, cell_type)
}

get_logFC_wide <- function(pathway_genes, datalist, ordered_genes) {
  bind_rows(lapply(names(datalist), function(ct) {
    datalist[[ct]] %>%
      filter(gene %in% pathway_genes) %>%
      select(gene, logFC, cell_type)
  })) %>%
    pivot_wider(names_from = cell_type, values_from = logFC, values_fill = 0) %>%
    mutate(gene = factor(gene, levels = ordered_genes)) %>%
    arrange(gene)
}

ordered_genes_list <- map(pathway_genes, get_ordered_genes, datalist, custom_order)
logFC_wide_list <- map2(pathway_genes, ordered_genes_list, ~ get_logFC_wide(.x, datalist, .y$gene))

# Add BBB full version
bbb_genes <- pathway_genes[["BBB"]]
bbb_present <- ordered_genes_list[["BBB"]]$gene
ordered_genes_list[["BBB_full"]] <- data.frame(gene = union(bbb_present, bbb_genes), cell_type = "all")

# Save results
saveRDS(ordered_genes_list, file = "../results/RNA/rds/ordered_genes_list.rds")
saveRDS(logFC_wide_list, file = "../results/RNA/rds/logFC_wide_list.rds")

logFC_wide <- lapply(names(logFC_wide_list), function(nm) {
  x <- logFC_wide_list[[nm]]
  x$process <- nm
  return(x)
})

write.csv(bind_rows(logFC_wide), "../results/RNA_ATAC_integration/tables/Genes_process_FC.csv", row.names = FALSE)
# Load links and filter peaks
links_df <- read.csv("../results/ATAC/data/endo_linkpeaks.csv")
links_to_join <- links_df %>%
  select(gene, peaks) %>%
  filter(peaks %in% peaks.v)

# Set display and file order for heatmaps
ordered_genes_list <- ordered_genes_list[c(11, 7, 6, 10, 4,  1, 13, 5, 12, 15, 14, 9, 8, 3, 2, 16)]
# Link DEGs to ATAC cCREs
deg_flox_linked <- map(ordered_genes_list, function(df) {
  df %>%
    left_join(log2cpm.scaled_rna, by = "gene") %>%
    left_join(links_to_join, by = "gene") %>%
    left_join(log2cpm.scaled_atac, by = "peaks") %>%
    drop_na()
})

# Save gene-peak tables for HOMER or export
dir.create("../results/RNA_ATAC_integration/pathways_gene_enhancer", recursive = TRUE, showWarnings = FALSE)
iwalk(deg_flox_linked, function(df, name) {
  write.csv(df %>% select(gene, peaks, cell_type),
            file = glue("../results/RNA_ATAC_integration/pathways_gene_enhancer/{name}_flox_gene_cCRE.csv"),
            row.names = FALSE)
})

# Save flat summary table of all gene-peak-process mappings
deg_flox_linked_noCts <- imap_dfr(deg_flox_linked, ~ {
  .x %>%
    select(-contains("Flox")) %>%
    mutate(process = .y)
})

# Rename cell types for display
celltype_levels <- c("Artery", "Capillary artery", "Capillary", "Capillary vein", "Large Vein")
celltype_map <- c("Artery" = "Artery", "CapArt" = "Capillary artery", "Cap" = "Capillary",
                  "CapVein" = "Capillary vein", "LargeVein" = "Large Vein")

deg_flox_linked_noCts_df <- deg_flox_linked_noCts %>%
  mutate(cell_names = recode(cell_type, !!!celltype_map),
         cell_names = factor(cell_names, levels = celltype_levels)) %>%
  select(cell_names, process, everything(), -cell_type)

# Summary tables
write.csv(deg_flox_linked_noCts_df, "../results/RNA_ATAC_integration/tables/process_genes_cCRE_names.csv", row.names = FALSE)

deg_summary <- deg_flox_linked_noCts_df %>%
  group_by(process) %>%
  summarise(n_genes = n_distinct(gene), n_peaks = n_distinct(peaks))

deg_summary_cell <- deg_flox_linked_noCts_df %>%
  group_by(process, cell_names) %>%
  summarise(n_genes = n_distinct(gene), n_peaks = n_distinct(peaks))

write.csv(deg_summary, "../results/RNA_ATAC_integration/tables/summary_process_genes_cCRE.csv", row.names = FALSE)
write.csv(deg_summary_cell, "../results/RNA_ATAC_integration/tables/summary_process_genes_cCRE_per_celltype.csv", row.names = FALSE)


# === Create output directory for heatmaps ===
heatmap_dir <- "../results/RNA_ATAC_integration/graphs/process_genes_enhancers_flox"
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

# === Annotation setup ===

EC_colors <- c("#247567", "#E2705B", "#6C7A99", "#B4698B", "#D8A21F")
names(EC_colors) <- c("Artery", "Capillary artery", "Capillary", "Capillary vein", "Large Vein")
my_colour <- list("Endothelial cell" = EC_colors)

top_annot_gene <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = "#C03830", col = NA)),
                                    `Endothelial cell` = rep(c("Artery", "Capillary artery", "Capillary", "Capillary vein", "Large Vein"), each = 3),
                                    col = my_colour,
                                    annotation_height = unit(c(2, 2), "mm"),  # Adjust heights for bars
                                    show_legend = c(FALSE, FALSE),
                                    show_annotation_name = FALSE,
                                    annotation_legend_param = list(title_gp = gpar(fontsize = 5, fontface = "bold"),   # Legend title appearance
                                                                   labels_gp = gpar(fontsize = 5, col = "black")))# 👈 Hide legend for Source, keep for cell types

top_annot_peak <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = "#6B6ACf", col = NA)),
                                    `Endothelial cell` = rep(c("Artery", "Capillary artery", "Capillary", "Capillary vein", "Large Vein"), each = 3),
                                    col = my_colour,
                                    annotation_height = unit(c(2, 2), "mm"),  # Adjust heights for bars
                                    show_legend = c(FALSE, FALSE),
                                    show_annotation_name = FALSE,
                                    annotation_legend_param = list(title_gp = gpar(fontsize = 5, fontface = "bold"),   # Legend title appearance
                                                                   labels_gp = gpar(fontsize = 5, col = "black")))# 👈 Hide legend for Source, keep for cell types

# === Heatmap titles for pathways ===
titles_heatmap <- c("Angiogenesis", "Transcription factors", "Inflammatory response", 
                    "Coagulation", "Cell metabolism", 
                   "Hypoxia", "Autophagy", "Transporters", 
                    "Chromatin remodeling", "Mitochondrion organization",  "Mitophagy", 
                    "Xenobiotic transport", "Glucose metabolism", "Flow sensing", "BBB", "BBB full")


names(titles_heatmap) <- names(deg_flox_linked)

# === Top marker genes per process for annotation ===

top1_per_cell_process <- bind_rows(logFC_wide) %>%
  filter(gene %in% deg_flox_linked_noCts$gene) %>% 
  # only keep genes in exactly one process
  group_by(gene) %>%
  filter(n_distinct(process) == 1) %>%
  ungroup() %>%
  # long form
  pivot_longer(cols      = Artery:LargeVein,
    names_to  = "cell_type",
    values_to = "logFC") %>%
  # first: pick the top row per (process × cell_type)
  group_by(process, cell_type) %>%
  slice_max(logFC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # then: for any gene that still appears more than once in a process,
  # keep only the one assignment where its logFC is highest
  group_by(process, gene) %>%
  slice_max(logFC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # finally, reorder & sort if you like
  mutate(process   = factor(process, levels = names(deg_flox_linked)),
    cell_type = factor(cell_type,
                       levels = c("Artery","CapArt","Cap","CapVein","LargeVein"))) %>%
  arrange(process, cell_type)

top1_per_cell_process

# === Plot heatmaps ===
letters_index <- letters[1:length(deg_flox_linked)] 
index <- 1
for (name in names(deg_flox_linked)[-16]) {
  gene_linkpeaks <- deg_flox_linked[[name]]
  full_name <- titles_heatmap[[name]]
  RNA_label_motif <- top1_per_cell_process %>% filter(process == name) %>% pull(gene)
  rownames_all <- gene_linkpeaks$gene
  gTF_positions <- match(RNA_label_motif, rownames_all)
  gTF_positions_anno <- gTF_positions[!is.na(gTF_positions)]
  
  ha_gTF <- rowAnnotation(motif_labels = anno_mark(at = gTF_positions_anno,
                                                   labels = RNA_label_motif[!is.na(gTF_positions)],
                                                   labels_gp = gpar(fontsize = 10, col = "black"),
                                                   padding = unit(1, "mm"),
                                                   extend = unit(0.2, "mm")))

  col_gene <- colorRamp2(seq(-1.9, 1.9, length.out = 100), plasma(100))
  
  hmap_gene <- Heatmap(gene_linkpeaks[ , c(3:17)],
                                        col = col_gene,
                                        cluster_columns = FALSE, cluster_rows = FALSE,
                                        column_title = if (index == 1) "Gene expression" else NULL,
                                        column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                                        top_annotation = top_annot_gene,
                                        name = "Z-score RNA\nlog2(CPM +1)",
                                        show_row_names = FALSE, show_column_names = FALSE,
                                        show_heatmap_legend = FALSE,
                                        border = "black", heatmap_width = unit(2, "in"),
                                        use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Hanning")
 

  col_peaks <- colorRamp2(seq(-1.6, 1.6, length.out = 100), viridis(100))
  
  hmap_peak <- Heatmap(gene_linkpeaks[ , c(19:33)],
                       col = col_peaks,
                       cluster_columns = FALSE, cluster_rows = FALSE,
                       right_annotation = ha_gTF,
                       column_title = if (index == 1) "Putative enhancers" else NULL,
                       column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                       top_annotation = top_annot_peak,
                       name = "Z-score ATAC\nlog2(CPM +1)",
                       show_heatmap_legend = FALSE,
                       border = "black", width = unit(1.8, "in"),
                       show_row_names = FALSE, show_column_names = FALSE,
                       use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Hanning")
  
  combined <- hmap_gene + hmap_peak
  matrix_length <- nrow(gene_linkpeaks)
  
  png(filename = glue("{heatmap_dir}/{letters_index[index]}_{name}.png"), 
      width = 5, height = 0.00000517 * matrix_length + 1.387875, units = "in", res = 500)
  draw(combined, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
       column_title = paste0(full_name, "         "), column_title_gp = gpar(fontsize = 13, fontface = "bold"))
  dev.off()
  
  index <- index + 1
}

EC_legend <- Legend(labels = c("Artery", "Capillary artery", "Capillary", "Capillary vein", "Large Vein"), 
                    title = "Endothelial cell", 
                    legend_gp = gpar(fill = EC_colors),
                    title_gp = gpar(fontsize = 9.5, fontface = "bold"),   # Legend title appearance
                    labels_gp = gpar(fontsize = 9, col = "black"),
                    title_gap = unit(1, "mm"))


color_fun1 <- circlize::colorRamp2(c(seq(-2, 2, length.out = 100)), plasma(100))
color_fun2 <- circlize::colorRamp2(c(seq(-2, 2, length.out = 100)), viridis(100))

rna_legend <- Legend(col_fun = color_fun1, 
                     title = "Z-score RNA\nlog2(CPM +1)", 
                     at = c(-2, -1, 0, 1, 2),  # Define ticks on the legend
                     labels = c("-2", "-1", "0", "1", "2"),
                     title_gp = gpar(fontsize = 8, fontface = "bold"),   # Legend title appearance
                     labels_gp = gpar(fontsize = 8, col = "black"),
                     legend_height = unit(2, "cm"),
                     title_gap = unit(1, "mm"))

atac_legend <- Legend(col_fun = color_fun2, 
                      title = "Z-score ATAC\nlog2(CPM +1)", 
                      at = c(-2, -1, 0, 1, 2),  # Define ticks on the legend
                      labels = c("-2", "-1", "0", "1", "2"),
                      title_gp = gpar(fontsize = 8, fontface = "bold"),   # Legend title appearance
                      labels_gp = gpar(fontsize = 8, col = "black"),
                      legend_height = unit(2, "cm"),
                      title_gap = unit(1, "mm"))

combined_legend <- packLegend(EC_legend, rna_legend, atac_legend, direction = "horizontal")
png(filename = glue("{heatmap_dir}/legend.png"), width = 5, height = 1.5, units = "in", res = 500)
ComplexHeatmap::draw(combined_legend)
dev.off()

image_list <- list.files(glue("{heatmap_dir}"), pattern = "^\\D_.*$", full.names = TRUE)
legend_path <- list.files(glue("{heatmap_dir}"), pattern = "legend", full.names = TRUE)


all_images <- magick::image_read(c(image_list, legend_path))
all_images_selected <- magick::image_read(c(image_list[1:9], legend_path))
all_images_selected_noflow <- magick::image_read(c(image_list[-14], legend_path))

# Combine all images vertically
combined_all_vertical <- magick::image_append(all_images, stack = TRUE)
combined_all_vertical_sel <- magick::image_append(all_images_selected, stack = TRUE)
combined_all_vertical_sel_noflow <- magick::image_append(all_images_selected_noflow, stack = TRUE)

magick::image_write(combined_all_vertical, path = glue("{heatmap_dir}/all_combined_vertical.png"))
magick::image_write(combined_all_vertical_sel, path = glue("{heatmap_dir}/all_combined_vertical_sel.png"))
magick::image_write(combined_all_vertical_sel_noflow, path = glue("{heatmap_dir}/all_combined_vertical_sel_no_flow.png"))


### files for Homer

enhancers_pathwayflox_folder <- "../results/RNA_ATAC_integration/Homer/pathway_enhancers"
enhancers_pathwayflox_celltype <- "../results/RNA_ATAC_integration/Homer/pathway_enhancers_celltype"


dir.create(enhancers_pathwayflox_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(enhancers_pathwayflox_celltype, recursive = TRUE, showWarnings = FALSE)

dirs_to_create <- paste0(enhancers_pathwayflox_celltype, "/", names(deg_flox_linked[-16]))
purrr::walk(dirs_to_create, ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))


peaks_process_celltype <- lapply(deg_flox_linked[-16], function(x){
  x %>% select(gene, peaks, cell_type)
})


for (loop_var in names(peaks_process_celltype)) {
  
  #getting the data
  tmp <- peaks_process_celltype[[loop_var]]
  if (nrow(tmp) == 0) next
  
  #getting peaks per pathway
  peaks_tmp <- tmp %>% 
    select("peaks") %>% 
    distinct(peaks, .keep_all = TRUE) %>% 
    mutate(to_split = peaks) %>% 
    separate(to_split, into = c("chrom", "chromStart", "chromEnd"), sep = "-") %>% 
    mutate(across(c(chromStart, chromEnd), as.integer)) %>%
    select(chrom, chromStart, chromEnd, peaks)
  
  #saving peaks per pathway for Homer
  write.table(x = peaks_tmp, file = glue("{enhancers_pathwayflox_folder}/{loop_var}_homer.txt"), sep = "\t", quote = F, row.names = F)
  
  #getting cell types cCRE-gene per pathways
  cell_types <- unique(tmp$cell_type)
  
  #process cell type specific cCRE
  for (cell in cell_types) {
    dup_count <- tmp %>% 
      filter(cell_type == cell) %>%
      summarize(duplicates = sum(duplicated(peaks))) %>%
      pull(duplicates)
    
    
    if (dup_count > 0) cat("⚠️ Duplicates found for", cell, "in", loop_var, ":", dup_count, "\n")
    
    #getting peaks per cell type in current pathway
    peaks_celltype_tmp <- tmp %>% 
      filter(cell_type == cell) %>% 
      distinct(peaks, .keep_all = TRUE) %>% 
      mutate(to_split = peaks) %>% 
      separate(to_split, into = c("chrom", "chromStart", "chromEnd"), sep = "-") %>% 
      mutate(across(c(chromStart, chromEnd), as.integer)) %>%
      select(chrom, chromStart, chromEnd, peaks)
    
    # Get unique chromosomes
    unique_chromosomes <- unique(peaks_celltype_tmp$chrom)
    
    # Logging duplicates removed
    if (dup_count > 0) {
      cat("⚡ Removed", dup_count, "duplicates for", cell, "in", loop_var, "\n")
    }
    cat("🧬 Unique chromosomes for", cell, "in", loop_var, ":", paste(unique_chromosomes, collapse = ", "), "\n\n")
    
    #saving for Homer
    write.table(x = peaks_celltype_tmp, file = glue("{enhancers_pathwayflox_celltype}/{loop_var}/{cell}_{loop_var}_homer.txt"), sep = "\t", quote = F, row.names = F)
    
  }
}


# Optionally save session info
sink("../results/logs/sessioninfo_11_RNA_ATAC_pathway_integration_Flox_BEC_subtypes.txt")
sessionInfo()
sink()

