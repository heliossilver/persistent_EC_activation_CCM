
pacman::p_load("tibble", "RColorBrewer", "pheatmap", "colorspace",
               "stringr","purrr", "viridis", "ggtext",
               "ComplexHeatmap", "circlize", "glue","patchwork",
               "dplyr", "ggplot2", "textshape", "paletteer", "tibble", "tidyr")

# Output folders----
tables_folder <- "../results/RNA_ATAC_integration/tables"
graphs_folder <- "../results/RNA_ATAC_integration/graphs"
enhancer_flox_dir <- "../results/RNA_ATAC_integration/Homer/enhancer_flox"
cCRE_KO_folder <- "../results/RNA_ATAC_integration/Homer/cCRE_KO"
enhancers_KO_folder <- "../results/RNA_ATAC_integration/Homer/enhancers_KO"

dir.create(tables_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(graphs_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(enhancer_flox_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cCRE_KO_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(enhancers_KO_folder, recursive = TRUE, showWarnings = FALSE)
#######Process Flox----

## Links----

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


links_df <- read.csv(file = "../results/ATAC/data/endo_linkpeaks.csv")


links_to_join <- links_df %>% 
  select(gene, peaks) 

### genes list ----
files <- list.files(path = "../results/RNA/edgeR/deg_difexp_tables/", pattern = "_vs_rest_diffexp.csv")
samples <- gsub("_vs_rest_diffexp.csv", "", files[c(1, 3, 2, 4, 5)])


up_genes.l <- list()
for (sample in samples) {
  
  data <- read.csv(file = glue("../results/RNA/edgeR/deg_difexp_tables/{sample}_vs_rest_diffexp.csv")) %>% 
    filter(diffexpressed == "Up") %>% 
    select(gene, starts_with("logFC")) %>% 
    rename_with(~ "logFC", starts_with(paste0("logFC_", sample, "_vs_rest"))) %>% 
    mutate(cell_type = paste0(sample))
  
  up_genes.l[[sample]] <- data
  
}

up_genes.df <- bind_rows(up_genes.l)
up_genes.df <- up_genes.df %>% 
  arrange(desc(logFC)) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  mutate(cell_type = factor(cell_type, levels = c("Artery", "CapArt", "Cap", "CapVein", "LargeVein"))) %>%
  arrange(cell_type) 

up_genes <- pull(up_genes.df, var = "gene")

### peaks lists----

files <- list.files(path = "../results/ATAC/edgeR/dar_diffacc_tables/", pattern = "_vs_rest_diffacc.csv")
samples <- gsub("_vs_rest_diffacc.csv", "", files[c(1, 3, 2, 4, 5)])

up_peaks.l <- list()
for (sample in samples) {
  
  data <- read.csv(file = glue("../results/ATAC/edgeR/dar_diffacc_tables/{sample}_vs_rest_diffacc.csv")) %>%
    filter(diffaccessible == "Up") %>%
    select(peaks, starts_with("logFC")) %>%
    rename_with(~ "logFC", starts_with(paste0("logFC_", sample, "_vs_rest"))) %>%
    mutate(cell_type = paste0(sample))
  
  up_peaks.l[[sample]] <- data
  
}


up_peaks.df <- bind_rows(up_peaks.l)
up_peaks.df <- up_peaks.df %>%
  arrange(desc(logFC)) %>%
  distinct(peaks, .keep_all = TRUE)

up_peaks.df$cell_type <- factor(up_peaks.df$cell_type,
                                levels = c("Artery", "CapArt", "Cap", "CapVein", "LargeVein"))

# Order the data frame by cell_type
up_peaks.df <- up_peaks.df %>%
  arrange(cell_type)

up_peaks <- pull(up_peaks.df, var = "peaks")


# ### Joining DEG with peak CPMs through linked peaks-genes ----

links_to_join_filt <- links_to_join %>%
  filter(peaks %in% up_peaks)


gene_linkpeaks <- data.frame(gene = up_genes)
gene_linkpeaks <- gene_linkpeaks %>%
  left_join(log2cpm.scaled_rna, by = "gene") %>%
  left_join(links_to_join_filt, by = "gene") %>%
  #left_join(links_to_join, by = "gene") %>%
  left_join(log2cpm.scaled_atac, by = "peaks") %>% 
  drop_na(peaks)

gene_linkpeaks[is.na(gene_linkpeaks)] <- 0

gene_celltype <- up_genes.df %>% select(gene, cell_type)

gene_link_table <- gene_linkpeaks %>% 
  left_join(gene_celltype, by = "gene") %>% 
  select(peaks, gene, cell_type)

gene_link_table_filt <- gene_link_table %>% mutate(peak_gene_pairs = paste(peaks, gene, sep = "_")) %>% distinct(peak_gene_pairs, .keep_all = TRUE)

write.csv(gene_link_table, file = glue("{tables_folder}/gene_enhancer_flox_per_cell.csv"), row.names = FALSE)

### saving enhancer peaks for homer analysis
# Total duplicate count across all cell types
total_duplicates <- sum(duplicated(gene_link_table$peaks))
cat("🔍 Total duplicated peaks across all cell types:", total_duplicates, "\n")

gene_link_table_unique <- gene_link_table %>%
  group_by(peaks) %>%
  filter(n() == 1) %>%  # Keep only peaks that appear once across all cell types
  ungroup()

# Check how many rows remain and how many were removed
cat("✅ Remaining unique peaks:", nrow(gene_link_table_unique), "\n")
cat("❌ Removed peaks shared across cell types:", nrow(gene_link_table) - nrow(gene_link_table_unique), "\n")

#check number per cell type
gene_link_table_unique %>% group_by(cell_type) %>% 
  summarise(number = n()) 

cell_types <- unique(gene_link_table$cell_type)

#process cell type specific cCRE
for (cell in cell_types) {
  
  #getting peaks per cell type in current pathway
  peaks_celltype_tmp <- gene_link_table_unique %>% 
    filter(cell_type == cell) %>% 
    mutate(to_split = peaks) %>% 
    separate(to_split, into = c("chrom", "chromStart", "chromEnd"), sep = "-") %>% 
    mutate(across(c(chromStart, chromEnd), as.integer)) %>%
    select(chrom, chromStart, chromEnd, peaks)
  # Progress message
  cat("📝 Processing", cell, "with", nrow(peaks_celltype_tmp), "unique peaks.\n")
  
  # Check if the data frame is empty before saving
  if (nrow(peaks_celltype_tmp) > 0) {
    write.table(x = peaks_celltype_tmp, 
                file = glue("{enhancer_flox_dir}/{cell}_markers_homer.txt"), 
                sep = "\t", quote = F, row.names = F)
    cat("✅ Saved:", glue("{enhancer_flox_dir}/{cell}_markers_homer.txt"), "\n\n")
  } else {
    cat("⚠️ No unique peaks for", cell, "- file not saved.\n\n")
  }
}

### process for heatmaps
gene_linkpeaks_split <- gene_linkpeaks %>% 
  left_join(gene_celltype, by = "gene") %>% 
  group_by(cell_type) %>% 
  summarise(number = n())  %>% 
  pull(number)

number_genes_bec <- gene_linkpeaks %>% 
  left_join(gene_celltype, by = "gene") %>% 
  select(gene, cell_type) %>%
  group_by(cell_type) %>%
  summarise(unique_genes = n_distinct(gene)) %>% 
  pull(unique_genes)

number_peaks_bec <- gene_linkpeaks %>% 
  left_join(gene_celltype, by = "gene") %>% 
  select(peaks, cell_type) %>%
  group_by(cell_type) %>%
  summarise(unique_peaks = n_distinct(peaks)) %>% 
  pull(unique_peaks)


### heatmap genes-cCRE subtype----
annot <- data.frame(`Endothelial cell` = rep(rep(c("Artery", "Capillary artery", "Capillary", "Capillary vein", "Large vein"),
                                                 c(3, 3, 3, 3, 3)), 2),
                    #source = rep(c("RNAseq", "ATACseq"), c(15,15)),
                    row.names = colnames(gene_linkpeaks[c(2:16, 18:32)]))
colnames(annot) <- "Endothelial cell"
annot$`Endothelial cell` <- factor(annot$`Endothelial cell`, levels = unique(annot$`Endothelial cell`))
#annot$Source <- factor(annot$Source, levels = unique(annot$Source))



# Generate the color palette
EC_colors <-c("#247567", 
              "#E2705B",
              "#6C7A99",
              "#B4698B",
              "#D8A21F")

names(EC_colors) <-  paste0(unique(annot$`Endothelial cell`))


my_colour <-  list("Endothelial cell"= EC_colors)


gaps_row_gene <- cumsum(gene_linkpeaks_split[1:4])

top_annot_gene <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = "#C03830", col = NA)),
                                    `Endothelial cell` = annot$`Endothelial cell`[1:15],
                                    col = list(`Endothelial cell` = EC_colors),
                                    annotation_height = unit(c(2.5, 2.5), "mm"),  # Adjust heights for bars
                                    show_legend = c(FALSE, TRUE),
                                    show_annotation_name = FALSE,
                                    annotation_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),   # Legend title appearance
                                                                   labels_gp = gpar(fontsize = 8, col = "black")))# 👈 Hide legend for Source, keep for cell types


left_annotation_gene <- rowAnnotation(Number_genes = anno_block(labels = number_genes_bec, labels_rot = 0,
                                                                labels_gp = gpar(col = "black", fontsize = 6, fontface = "bold")),
                                      width = unit(4.4, "mm"), height = unit(3.2, "mm"),
                                      show_annotation_name = TRUE)

hmap_flox_gene <- ComplexHeatmap::pheatmap(gene_linkpeaks[ , c(2:16)],  
                                           color = plasma(100),
                                           cluster_cols = F, 
                                           cluster_rows = F, 
                                           breaks = seq(-2, 2, length.out = 100),
                                           main = "Gene\nexpression", border_color = NA,
                                           top_annotation = top_annot_gene,
                                           name = "Z-score RNA\nlog2(CPM +1)", 
                                           heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),   # Legend title appearance
                                                                       labels_gp = gpar(fontsize = 6, col = "black"),
                                                                       legend_height = unit(1.5, "cm")), 
                                           show_colnames = F,
                                           show_rownames = F, 
                                           fontsize = 8,
                                           left_annotation = left_annotation_gene,
                                           border = "black",
                                           annotation_colors = my_colour, 
                                           gaps_row = gaps_row_gene, 
                                           use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Hanning")



top_annot_peak <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = "#6B6ACf", col = NA)),  # 💡 col = NA removes the border
                                    `Endothelial cell` = annot$`Endothelial cell`[16:30],
                                    col = list(`Endothelial cell` = EC_colors),
                                    annotation_height = unit(c(2.5, 2.5), "mm"),
                                    show_legend = c(FALSE, TRUE),
                                    show_annotation_name = FALSE,
                                    annotation_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),   # Legend title appearance
                                                                   labels_gp = gpar(fontsize = 8, col = "black")))  # ✅ Hides annotation name

left_annotation_peak <- rowAnnotation(Number_peaks = anno_block(labels = number_peaks_bec, labels_rot = 0,
                                                                labels_gp = gpar(col = "black", fontsize = 6, fontface = "bold")),
                                      width = unit(4.7, "mm"), height = unit(3.2, "mm"),
                                      show_annotation_name = TRUE)

hmap_flox_peak <- ComplexHeatmap::pheatmap(gene_linkpeaks[ , c(18:32)],
                                           col = viridis(100),
                                           cluster_cols = F,
                                           cluster_rows = F,
                                           breaks = seq(-1.5, 1.5, length.out = 100),
                                           main = "Putative\nenhancers", border_color = NA,
                                           name = "Z-score ATAC\nlog2(CPM +1)", annotation_names_col = F,
                                           show_colnames = F, 
                                           heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),   # Legend title appearance
                                                                       labels_gp = gpar(fontsize = 6, col = "black"),
                                                                       legend_height = unit(1.5, "cm")), 
                                           fontsize = 8, 
                                           show_rownames = F, 
                                           top_annotation = top_annot_peak,
                                           #right_annotation = width_annotation,
                                           left_annotation = left_annotation_peak,
                                           annotation_colors = my_colour, 
                                           gaps_row = gaps_row_gene, 
                                           border = "black",
                                           use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Hanning")

combined_heatmaps <- hmap_flox_gene + hmap_flox_peak

# Draw them together

png(filename = glue("{graphs_folder}/flox_cCRE_gene_2hmaps.png"), width = 3.5, height = 4, units = "in", res = 300)
ComplexHeatmap::draw(combined_heatmaps, heatmap_legend_side = "left",
                     annotation_legend_side = "left", merge_legends = TRUE)
dev.off()


### Genomic annotation of flox accesible enhancers
peaks_present <- gene_linkpeaks %>% pull(peaks) %>% unique()

annot_peaks <- read_tsv("../results/ATAC/data/annotated_union_clean_sorted.peakset.txt") %>% 
  rename(peaks = starts_with("PeakID")) %>% 
  select(peaks, Annotation) %>% 
  mutate(Annotation_final = case_when(str_detect(Annotation, regex("promoter", ignore_case = TRUE)) ~ "Promoter-TSS",
                                      str_detect(Annotation, regex("5' UTR", ignore_case = TRUE)) ~ "5' UTR",
                                      str_detect(Annotation, regex("3' UTR", ignore_case = TRUE)) ~ "3' UTR",
                                      str_detect(Annotation, regex("exon", ignore_case = TRUE)) ~ "Exon",
                                      str_detect(Annotation, regex("Intron", ignore_case = TRUE)) ~ "Intron",
                                      str_detect(Annotation, regex("intergenic", ignore_case = TRUE)) ~ "Distal\nIntergenic",
                                      str_detect(Annotation, regex("TTS", ignore_case = TRUE)) ~ "TTS",
                                      TRUE ~ Annotation))

                                      
annotation_counts <- annot_peaks %>% 
  filter(peaks %in% peaks_present) %>% 
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
  theme_void() + ggtitle("Genomic annotation of putative enhancers <br><i>PDCD10</i><sup>fl/fl</sup> ") +
  theme(plot.title = element_markdown(hjust = 0.5, face = "bold"), 
        legend.title = element_blank()) + annotate("text", x = 0.95, y = 0, label = "", size = 1)

ggsave(filename = glue("{graphs_folder}/putative_enhancer_annotation_flox.svg"), plot = annot_peaks_flox_enhancers, 
       device = "svg", width = 4.8, height = 2.5, units = "in", bg = "white")

levels(gene_link_table$cell_type)

WT_gene_info_list <- list()
WT_peak_info_list <- list()
link_table_list <- list()
for (loop_var in levels(gene_link_table$cell_type)) {
  
  tmp_data <- gene_link_table %>% 
    filter(cell_type == loop_var)
  
  tmp_rna <- read.csv(glue("../results/RNA/edgeR/deg_difexp_tables/{loop_var}_vs_rest_diffexp_full.csv")) %>% 
    select(gene, starts_with("logFC"), starts_with("FDR"), diffexpressed) %>%
    rename_with(~glue("logFC_{loop_var}_vs_rest_gene"), starts_with("logFC")) %>% 
    rename_with(~glue("FDR_{loop_var}_vs_rest_gene"), starts_with("FDR")) %>% 
    rename_with(~glue("diffexpressed_{loop_var}"), starts_with("diff"))
  
  tmp_atac <- read.csv(glue("../results/ATAC/edgeR/dar_diffacc_tables/{loop_var}_KO_vs_Flox_diffacc_full.csv")) %>% 
    select(peaks, starts_with("logFC"), starts_with("FDR"), diffaccessible) %>% 
    rename_with(~glue("logFC_{loop_var}_vs_rest_peak"), starts_with("logFC")) %>% 
    rename_with(~glue("FDR_{loop_var}_vs_rest_peak"), starts_with("FDR")) %>% 
    rename_with(~glue("diffaccessible_{loop_var}"), starts_with("diff"))
  
  
  WT_gene_info_list[[loop_var]] <- tmp_rna
  WT_peak_info_list[[loop_var]] <- tmp_atac
  link_table_list[[loop_var]] <- tmp_data
}


all_rna_info <- purrr::reduce(WT_gene_info_list, full_join, by = "gene") %>%
  mutate(across(starts_with("logFC"), ~ replace_na(.x, 0)),
         across(starts_with("FDR"), ~ replace_na(.x, 1)),
         across(starts_with("diff"), ~ replace_na(.x, "No")))

all_atac_info <- purrr::reduce(WT_peak_info_list, full_join, by = "peaks") %>%
  mutate(across(starts_with("logFC"), ~ replace_na(.x, 0)),
         across(starts_with("FDR"), ~ replace_na(.x, 1)),
         across(starts_with("diff"), ~ replace_na(.x, "No")))

all_links <- bind_rows(link_table_list) %>% distinct()

full_WT_cCRE_FC_FDR <- all_links %>%
  left_join(all_rna_info, by = "gene") %>%
  left_join(all_atac_info, by = "peaks") %>% select(-cell_type)

cell_order <- c("Artery", "CapArt", "Cap", "CapVein", "LargeVein")

meta_cols <- c("gene", "peaks")

cell_cols <- setdiff(colnames(full_WT_cCRE_FC_FDR), meta_cols)

cell_cols_ordered <- unlist(
  lapply(cell_order, function(ct) {
    cell_cols[str_detect(cell_cols, paste0("_", ct, "($|_)"))]
  })
)

full_WT_cCRE_FC_FDR <- full_WT_cCRE_FC_FDR %>%
  select(all_of(meta_cols), all_of(cell_cols_ordered))



write.csv(x = full_WT_cCRE_FC_FDR, file = glue("{tables_folder}/gene_enhancer_flox_per_cell_FC.csv"), row.names = FALSE)



#######Process KO----


## Links----


links_df <- read.csv(file = "../results/ATAC/data/endo_linkpeaks.csv")
links_to_join <- links_df %>% 
  select(gene, peaks) 


files_rna <- gsub("_vs_Flox_diffexp.csv", "", list.files("../results/RNA/edgeR/deg_difexp_tables/", pattern = "_vs_Flox_diffexp.csv"))
files_atac <- gsub("_vs_Flox_diffacc.csv", "", list.files("../results/ATAC/edgeR/dar_diffacc_tables/", pattern = "_vs_Flox_diffacc.csv"))

samples <- if (all(files_rna %in% files_atac)) {
  unique(files_rna[c(1, 3, 2, 4, 5)])
} else {
  stop("Mismatch between RNA and ATAC samples.")
}

cell_types <- gsub("_KO", "", samples)

#atac
log2cpm_atac <- read.csv(file = "../results/ATAC/data/log2cpm.csv")


#rna
log2cpm_rna <- read.csv(file = "../results/RNA/data/log2cpm.csv", check.names = F) %>% column_to_rownames("gene") 


# Generate the color palette
EC_colors <-c("#247567", 
              "#E2705B",
              "#6C7A99",
              "#B4698B",
              "#D8A21F")

names(EC_colors) <-  cell_types

geno_colors <- c("#4DBBD5", "#F5AE97")
names(geno_colors) <- c("Flox", "KO")

my_colour <- list(`Endothelial cell` = EC_colors,
                  Genotype = geno_colors)

new_names <- c("Artery", "Capillary artery", "Capillary", "Capillary vein", "Large vein")
names(new_names) <- cell_types


col_names <- gsub("_Flox\\d|_KO\\d|_ATAC|_RNA","",colnames(log2cpm_atac[-1])) %>% unique()
names(col_names) <- cell_types

cts_rna <- list()
cts_atac <- list()
cts_linked <- list()
for (cell_type in cell_types) {
  
  tmp_rna <- read.csv(glue("../results/RNA/edgeR/deg_difexp_tables/{cell_type}_KO_vs_Flox_diffexp.csv")) %>% 
    select(gene, diffexpressed)
  tmp_atac <- read.csv(glue("../results/ATAC/edgeR/dar_diffacc_tables/{cell_type}_KO_vs_Flox_diffacc.csv")) %>% 
    filter(peaks %in% links_to_join$peaks) %>% select(peaks, diffaccessible)
  
  links_to_join_filt <-links_to_join %>% 
    filter(peaks %in% tmp_atac$peaks)
  
  col_name <- col_names[[cell_type]]
  
  selected_col_r <- colnames(log2cpm_rna)[str_detect(colnames(log2cpm_rna), paste0("^", col_name, "_[A-Z]"))]
  
  cts_rna_tmp <- log2cpm_rna %>% 
    select(all_of(selected_col_r)) %>% 
    rownames_to_column("gene")
  
  selected_col_a <- colnames(log2cpm_atac)[str_detect(colnames(log2cpm_atac), paste0("^", col_name))]
  
  cts_atac_tmp <- log2cpm_atac %>% 
    filter(peaks %in% tmp_atac$peaks) %>% 
    select(peaks, all_of(selected_col_a)) 
  
  cts_final <- tmp_rna %>% 
    left_join(cts_rna_tmp, by = "gene") %>% 
    left_join(links_to_join_filt, by = "gene") %>% 
    left_join(cts_atac_tmp, by = "peaks") %>% 
    drop_na() %>% select(-diffexpressed)
  
  cts_rna_link <- cts_final %>% 
    select(gene) %>% left_join(tmp_rna, by = "gene") %>% select(gene, diffexpressed) %>% 
    mutate(cell_type = cell_type) %>% distinct(gene, .keep_all = TRUE)
  
  cts_atac_link <- cts_final %>% 
    select(peaks) %>% left_join(tmp_atac, by = "peaks") %>% select(peaks, diffaccessible) %>% 
    mutate(cell_type = cell_type) %>% distinct(peaks, .keep_all = TRUE)
  
  cts_rna[[cell_type]] <- cts_rna_link
  cts_atac[[cell_type]] <- cts_atac_link
  cts_linked[[cell_type]] <- cts_final
  # 
  rm(tmp_rna, tmp_atac, cts_rna_link, cts_atac_link, cts_rna_tmp, cts_atac_tmp)
  # 
  ec_color <- EC_colors[[cell_type]]
  # 
  top_annot_gene <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = "#C03830", col = NA)),
                                      Genotype = rep(c("Flox", "KO"), c(3, 3)),
                                      `Endothelial cell` = rep(cell_type, 6),
                                      col = my_colour,
                                      annotation_height = unit(c(2.5, 2.5, 2.5), "mm"),  # Adjust heights for bars
                                      show_legend = c(FALSE, TRUE, FALSE),
                                      show_annotation_name = FALSE,
                                      annotation_legend_param = list(title_gp = gpar(fontsize = 5.5, fontface = "bold"),   # Legend title appearance
                                                                     labels_gp = gpar(fontsize = 6, col = "black")))# 👈 Hide legend for Source, keep for cell types
  
  
  hmap_flox_gene <- ComplexHeatmap::pheatmap(cts_final[ , c(2:7)],
                                             color = plasma(100), scale = "row",
                                             cluster_cols = F,
                                             cluster_rows = F,
                                             breaks = seq(-2, 2, length.out = 100),
                                             main = "Gene\nexpression", border_color = NA,
                                             top_annotation = top_annot_gene,
                                             name = "Z-score RNA\nlog2(CPM +1)",
                                             heatmap_legend_param = list(title_gp = gpar(fontsize = 5, fontface = "bold"),   # Legend title appearance
                                                                         labels_gp = gpar(fontsize = 6, col = "black"),
                                                                         legend_height = unit(1, "cm")),
                                             show_colnames = F,
                                             show_rownames = F,
                                             fontsize = 7,
                                             border = "black",
                                             annotation_colors = my_colour,
                                             use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Hanning")
  
  
  
  top_annot_peak <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = "#6B6ACf", col = NA)),  # 💡 col = NA removes the border
                                      Genotype = rep(c("Flox", "KO"), c(3, 3)),
                                      `Endothelial cell` = rep(cell_type, 6),
                                      col = my_colour,
                                      annotation_height = unit(c(2.5, 2.5, 2.5), "mm"),  # Adjust heights for bars
                                      show_legend = c(FALSE, TRUE, FALSE),
                                      show_annotation_name = FALSE,
                                      annotation_legend_param = list(title_gp = gpar(fontsize = 5, fontface = "bold"),   # Legend title appearance
                                                                     labels_gp = gpar(fontsize = 6, col = "black")))# 👈 Hide legend for Source, keep for cell types
  
  hmap_flox_peak <- ComplexHeatmap::pheatmap(cts_final[ , c(9:14)],
                                             col = viridis(100),
                                             cluster_cols = F, scale = "row",
                                             cluster_rows = F,
                                             breaks = seq(-2, 2, length.out = 100),
                                             main = "Putative\nenhancers", border_color = NA,
                                             name = "Z-score ATAC\nlog2(CPM +1)", annotation_names_col = F,
                                             show_colnames = F,
                                             heatmap_legend_param = list(title_gp = gpar(fontsize = 5, fontface = "bold"),   # Legend title appearance
                                                                         labels_gp = gpar(fontsize = 6, col = "black"),
                                                                         legend_height = unit(1, "cm")),
                                             fontsize = 7,
                                             show_rownames = F,
                                             top_annotation = top_annot_peak,
                                             annotation_colors = my_colour,
                                             border = "black",
                                             use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Hanning")
  #hmap_flox_peak
  
  combined_heatmaps <- hmap_flox_gene + hmap_flox_peak
  new_name <- new_names[[cell_type]]
  
  # Draw them together
  
  png(filename = glue("../results/RNA_ATAC_integration/graphs/{cell_type}_KO_cCRE_gene_2hmaps.png"), width = 1.5, height = 3.5, units = "in", res = 300)

  ComplexHeatmap::draw(combined_heatmaps, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
                       column_title_gp = gpar(fontsize = 12, fontface = "bold", col = ec_color),
                       column_title = glue("{new_name}"))
  dev.off()
  
  rm(cts_final, top_annot_gene, top_annot_peak, hmap_flox_gene, hmap_flox_peak)
  
}


EC_colors_lgnd <- EC_colors
names(EC_colors_lgnd) <- new_names

EC_legend <- Legend(labels = new_names, 
                    title = "Endothelial cell", 
                    legend_gp = gpar(fill = EC_colors),
                    title_gp = gpar(fontsize = 7, fontface = "bold"),   # Legend title appearance
                    labels_gp = gpar(fontsize = 7, col = "black"),
                    title_gap = unit(1, "mm"))

genotype_legend <- Legend(labels = names(geno_colors), 
                          title = "Genotype", 
                          legend_gp = gpar(fill = geno_colors),
                          title_gp = gpar(fontsize = 7, fontface = "bold"),   # Legend title appearance
                          labels_gp = gpar(fontsize = 7, col = "black"),
                          title_gap = unit(1, "mm"))


color_fun1 <- circlize::colorRamp2(c(seq(-2, 2, length.out = 100)), plasma(100))
color_fun2 <- circlize::colorRamp2(c(seq(-2, 2, length.out = 100)), viridis(100))

rna_legend <- Legend(col_fun = color_fun1, 
                     title = "Z-score RNA\nlog2(CPM +1)", 
                     at = c(-2,  0, 2),  # Define ticks on the legend
                     labels = c("-2",  "0", "2"),
                     title_gp = gpar(fontsize = 7, fontface = "bold"),   # Legend title appearance
                     labels_gp = gpar(fontsize = 7, col = "black"),
                     legend_height = unit(1, "cm"),
                     title_gap = unit(1, "mm"))

atac_legend <- Legend(col_fun = color_fun2, 
                      title = "Z-score ATAC\nlog2(CPM +1)", 
                      at = c(-2,  0,  2),  # Define ticks on the legend
                      labels = c("-2",  "0", "2"),
                      title_gp = gpar(fontsize = 7, fontface = "bold"),   # Legend title appearance
                      labels_gp = gpar(fontsize = 7, col = "black"),
                      legend_height = unit(1, "cm"),
                      title_gap = unit(1, "mm"))

combined_legend <- packLegend(EC_legend, genotype_legend, rna_legend, atac_legend, direction = "vertical")
png(filename = glue("{graphs_folder}/legend.png"), width = 1, height = 3.5, units = "in", res = 300)
ComplexHeatmap::draw(combined_legend)
dev.off()


#combined image ----
# Read images back and combine them horizontally


image_list <- list.files(glue("{graphs_folder}/"), pattern = "legend|2hmaps\\.png$", full.names = TRUE)

image_list_ordered <- image_list[c(7, 1, 3, 2, 4, 6)]

all_images <- magick::image_read(image_list_ordered)

# Combine all images vertically
combined_all <- magick::image_append(all_images)

magick::image_write(combined_all, path = glue("{graphs_folder}/all_combined.png"))




### saving for homer and full table----

gene_cCRE_KO <- lapply(cts_linked, function(x){
  x %>% select(gene, peaks)
}
)

gene_cCRE_KO_name <- Map(function(df, name) {
  df %>% 
    select(gene, peaks) %>%
    mutate(cell_type = name)
}, cts_linked, names(cts_linked))

write.csv(x = bind_rows(gene_cCRE_KO_name), file = glue("{tables_folder}/Enhancers_KO_per_cell.csv"), row.names = FALSE)
# gene_cCRE_KO <- read.csv(glue("{tables_folder}/Enhancers_KO_per_cell.csv")) %>%
#   mutate(cell_type = factor(cell_type, levels = c("Artery", "CapArt", "Cap", "CapVein", "LargeVein"))) %>%
#   split(.$cell_type) %>%
#   map(~select(.x, -cell_type))


KO_gene_info_list <- list()
KO_peak_info_list <- list()
link_table_list <- list()
for (loop_var in names(gene_cCRE_KO)) {
  
  tmp_data <- gene_cCRE_KO[[loop_var]]
  
  tmp_rna <- read.csv(glue("../results/RNA/edgeR/deg_difexp_tables/{loop_var}_KO_vs_Flox_diffexp_full.csv")) %>% 
    select(gene, starts_with("logFC"), starts_with("FDR"), diffexpressed) %>%
    rename_with(~glue("logFC_{loop_var}_KOvsFl_gene"), starts_with("logFC")) %>% 
    rename_with(~glue("FDR_{loop_var}_KOvsFl_gene"), starts_with("FDR")) %>% 
    rename_with(~glue("diffexpressed_{loop_var}_KO"), starts_with("diff"))
  
  tmp_atac <- read.csv(glue("../results/ATAC/edgeR/dar_diffacc_tables/{loop_var}_KO_vs_Flox_diffacc_full.csv")) %>% 
    select(peaks, starts_with("logFC"), starts_with("FDR"), diffaccessible) %>% 
    rename_with(~glue("logFC_{loop_var}_KOvsFl_peak"), starts_with("logFC")) %>% 
    rename_with(~glue("FDR_{loop_var}_KOvsFl_peak"), starts_with("FDR")) %>% 
    rename_with(~glue("diffaccessible_{loop_var}_KO"), starts_with("diff"))
  
  data_homer <- tmp_atac %>% 
    select(peaks, starts_with("diffaccessible")) %>% 
    filter(.[ ,2] == "Up") %>% 
    separate(peaks, into = c("chrom", "chromStart", "chromEnd"), sep = "-", remove = FALSE) %>% 
    select(chrom, chromStart, chromEnd, peaks)
  
  data_homer_enhancers <- filter(data_homer, peaks %in% tmp_data$peaks)
  # 
  write.table(x = data_homer, file = glue("{cCRE_KO_folder}/{loop_var}_KO_homer.txt"),
              sep = "\t", quote = F, row.names = F)
  write.table(x = data_homer_enhancers, file = glue("{enhancers_KO_folder}/{loop_var}_KO_homer.txt"),
              sep = "\t", quote = F, row.names = F)

  KO_gene_info_list[[loop_var]] <- tmp_rna
  KO_peak_info_list[[loop_var]] <- tmp_atac
  link_table_list[[loop_var]] <- tmp_data
}


all_rna_info <- purrr::reduce(KO_gene_info_list, full_join, by = "gene") %>%
  mutate(across(starts_with("logFC"), ~ replace_na(.x, 0)),
         across(starts_with("FDR"), ~ replace_na(.x, 1)),
         across(starts_with("diff"), ~ replace_na(.x, "No")))

all_atac_info <- purrr::reduce(KO_peak_info_list, full_join, by = "peaks") %>%
  mutate(across(starts_with("logFC"), ~ replace_na(.x, 0)),
         across(starts_with("FDR"), ~ replace_na(.x, 1)),
         across(starts_with("diff"), ~ replace_na(.x, "No")))

all_links <- bind_rows(link_table_list) %>% distinct()

full_KO_cCRE_FC_FDR <- all_links %>%
  left_join(all_rna_info, by = "gene") %>%
  left_join(all_atac_info, by = "peaks")

write.csv(x = full_KO_cCRE_FC_FDR, file = glue("{tables_folder}/Enhancers_KO_per_cell_FC.csv"), row.names = FALSE)

full_KO_cCRE_FC_FDR <- read.csv(glue("{tables_folder}/Enhancers_KO_per_cell_FC.csv"))

full_KO_cCRE_FC_strict <- full_KO_cCRE_FC_FDR %>%
  filter(if_all(matches("^FDR_.*_gene$"), ~ .x <= 0.05) &
           if_all(matches("^FDR_.*_peak$"), ~ .x <= 0.05))

write.csv(x = full_KO_cCRE_FC_strict, file = glue("{tables_folder}/enhancers_KO_per_cell_in_all_cells.csv"), row.names = FALSE)


## Genomic annotation of accesible enhancer KO
all_peaks_enhancersKO <- map(cts_atac, ~ .x$peaks) %>% unlist() %>% unique()
annotation_counts <- annot_peaks %>% 
  filter(peaks %in% all_peaks_enhancersKO) %>% 
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

annot_peaks_KO_enhancers <- ggplot(annotation_counts, aes(x = "", y = Percentage, fill = label)) +
  geom_bar(width = 1, stat = "identity", color = "white" ) +
  coord_polar("y") + 
  scale_fill_manual(values = color_mapping) +
  theme_void() + ggtitle("Genomic annotation of putative enhancers<br><i>PDCD10</i><sup>BECKO</sup> ") +
  theme(plot.title = element_markdown(hjust = 0.5, face = "bold"), 
        legend.title = element_blank()) + annotate("text", x = 0.95, y = 0, label = "", size = 1)

ggsave(filename = glue("{graphs_folder}/putative_enhancer_annotation_KO.svg"), plot = annot_peaks_KO_enhancers, 
       device = "svg", width = 4.8, height = 2.5, units = "in", bg = "white")

# ==== Reproducibility log ====
sink("../logs/sessioninfo_10_gene_peak_linkage_and_homer_files.txt")
sessionInfo()
sink()
