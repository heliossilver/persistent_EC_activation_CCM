
pacman::p_load("dplyr", "tibble", "tidyr", "ggplot2", "stringr", "purrr",
               "rtracklayer", "magick", "GenomicRanges", "readr",
               "vroom", "stringr", "viridis", "RColorBrewer", "glue", "EnrichedHeatmap", "forcats")

color_heatmaps <- brewer.pal(n = 8, name = "Set2")[c(1, 3, 4, 8)]
# "#66C2A5" "#8DA0CB" "#E78AC3" "#B3B3B3"

# 🔧 Function to normalize bigWig signals to matrices
normalize_bw_list <- function(bw_list, target, bin_val = 10, extend = 1000, mean_mode = "w0") {
  map2(
    bw_list,
    names(bw_list),
    ~ {
      cat("🔹 Normalizing:", .y, "\n")
      normalizeToMatrix(
        signal = .x,
        target = target,
        value_column = "score",
        extend = extend,
        mean_mode = mean_mode,
        w = bin_val
      )
    }
  ) %>%
    { set_names(., names(bw_list)) }  # ensure names are preserved
}
compute_quantiles <- function(mtx) {
  quantile(mtx, c(0.97, 0.99), na.rm = TRUE)
}

### Normalize To Matrix ----
junb_dep <- read.table("../results/Human_cells/ChipSeq/data/junb_positive_atac_up_peaks_less_than_500.bed", sep = "\t",
                       col.names =  c("chr", "start", "end")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

junb_indep <- read.table("../results/Human_cells/ChipSeq/data/junb_negative_atac_up_peaks_less_than_500.bed", sep = "\t",
                       col.names =  c("chr", "start", "end")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)


bw_file <- list.files("../results/Human_cells/D3_ATAC/bw/", full.names = T)
key_bw_names <- read.csv("../results/Human_cells/D3_ATAC/data/key_atac.csv", header = FALSE) 

bw_true_names <- setNames(key_bw_names$V1, key_bw_names$V2)

bw_names <- str_replace_all(basename(bw_file), c(bw_true_names, "_sorted_rmdup.bw" = ""))

# names(bw_file) <- bw_names

bw_list <- setNames(lapply(bw_file[c(1:3, 4:6, 29:30, 25, 33:34, 27, 7:9, 10:12, 31:32, 26, 35:36, 28)], import), 
                    bw_names[c(1:3, 4:6, 29:30, 25, 33:34, 27, 7:9, 10:12, 31:32, 26, 35:36, 28)]) #, 37, 39, 41, 43)])

dir.create("../results/Human_cells/bw_norm_mtx/junb_dep", showWarnings = FALSE, recursive = TRUE)
dir.create("../results/Human_cells/bw_norm_mtx/junb_indep", showWarnings = FALSE, recursive = TRUE)

full_junb <- c(junb_indep, junb_dep)


partition <- fct_inorder(c(rep(paste("JuNB-independent\n", length(junb_indep)), length(junb_indep)),
                          rep(paste("JUNB-dependent\n", length(junb_dep)), length(junb_dep))))


##JUNB chip ----
bw_file_chip <- list.files("../results/Human_cells/ChipSeq/bw/", full.names = T)[c(12, 11, 10, 9,
                                                                                   4, 3, 2, 1,
                                                                                   8, 7, 6, 5)]
bw_file_chip


bw_names_chip <- str_replace_all(basename(bw_file_chip), c("_sorted_rmdup.bw" = ""))

# names(bw_file) <- bw_names

bw_chip_list <- setNames(lapply(bw_file_chip, import), 
                         bw_names_chip)

# Normalize and save with correct filenames
norm_junbchip_mtx_l <- normalize_bw_list(bw_chip_list, target = full_junb)
# saveRDS(norm_junbchip_mtx_l, "../results/Human_cells/bw_norm_mtx/nnorm_junbchip_mtx_l.rds")
# norm_junbchip_mtx_l <- readRDS("../results/Human_cells/bw_norm_mtx/nnorm_junbchip_mtx_l.rds")

condition_names_chip <- c("PDCD10 WT\nVeh (36h)", "PDCD10 WT\nCCMenv (36h)", 
                          "PDCD10 KD\nVeh (36h)", "PDCD10 KD\nCCMenv (36h)")

titles_hmaps_chip <- setNames(rep(condition_names_chip, 3),
                              names(norm_junbchip_mtx_l))


mtx_junb_junb_wt_l <- map2(norm_junbchip_mtx_l[1:4], names(norm_junbchip_mtx_l)[1:4], 
                      ~{
                        .x %>%
                          as.data.frame() %>%
                          as.matrix()
                      }
)

mtx_h3ac_junb_wt_l <- map2(norm_junbchip_mtx_l[5:8], names(norm_junbchip_mtx_l)[5:8], 
                      ~{
                        .x %>%
                          as.data.frame() %>%
                          as.matrix()
                      }
)

mtx_h3me_junb_wt_l <- map2(norm_junbchip_mtx_l[9:12], names(norm_junbchip_mtx_l)[9:12], 
                      ~{
                        .x %>%
                          as.data.frame() %>%
                          as.matrix()
                      }
)

color_pal <- RColorBrewer::brewer.pal(12, "Paired")
color_pal_map <- setNames(c(color_pal[1:2], color_pal[5:6]),
                          condition_names_chip)

annot_junb_junb_l <- map2(titles_hmaps_chip, names(titles_hmaps_chip), 
                     ~{
                       HeatmapAnnotation(Source = anno_block(gp = gpar(fill = color_pal_map[[.x]], 
                                                                       col = "black")), 
                                         height = unit(2.5, "mm"),
                                         show_legend = FALSE)
                     }
)


### ordered heatmap----

row_order_ht <- row_order(EnrichedHeatmap(norm_junbchip_mtx_l[["junb_KD_CCM"]], pos_line = FALSE, show_heatmap_legend = FALSE,
                                          axis_name = "", split = partition,
                                          top_annotation = NULL))

order_for_heatmaps <- as.vector(unlist(row_order_ht))

# Normalize and save with correct filenames----
norm_mtx_junb_l <- normalize_bw_list(bw_list, full_junb)
#saveRDS(norm_mtx_junb_l, "../results/Human_cells/bw_norm_mtx/norm_mtx_junbdepind_l.rds")
#norm_mtx_junb_l <- readRDS("../results/Human_cells/bw_norm_mtx/norm_mtx_junbdepind_l.rds")

# Compute quantiles for norm----
quantiles_atac_junb <- lapply(norm_mtx_junb_l, compute_quantiles)
quantiles_junb_junb_dep <- lapply(norm_junbchip_mtx_l[1:4], compute_quantiles)
quantiles_h3ac_junb_dep <- lapply(norm_junbchip_mtx_l[5:8], compute_quantiles)
quantiles_h3me_junb_dep <- lapply(norm_junbchip_mtx_l[9:12], compute_quantiles)


condition_names <- rep(c("PDCD10 WT\nVeh (1h)", "PDCD10 WT\nCCMenv (1h)", 
                         "PDCD10 WT\nVeh (36h)", "PDCD10 WT\nCCMenv (36h)",
                         "PDCD10 KD\nVeh (1h)", "PDCD10 KD\nCCMenv (1h)", 
                         "PDCD10 KD\nVeh (36h)", "PDCD10 KD\nCCMenv (36h)"), each = 3)

titles_hmaps_junb <- setNames(condition_names,
                            names(norm_mtx_junb_l))




mtx_atac_junb_l <- map2(norm_mtx_junb_l, names(norm_mtx_junb_l), 
                      ~{
                        .x %>%
                          as.data.frame() %>%
                          as.matrix()
                      }
)

color_pal_atac_map <- setNames(c(rep(color_pal[1:2], each = 3), rep(color_pal[1:2], each = 3), 
                            rep(color_pal[5:6], each = 3), rep(color_pal[5:6], each = 3)),
                          condition_names)


annot_junb_l <- map2(titles_hmaps_junb, names(titles_hmaps_junb), 
                ~{
                  HeatmapAnnotation(Source = anno_block(gp = gpar(fill = color_pal_atac_map[[.x]], 
                                                                  col = "black")), 
                                    height = unit(2.5, "mm"),
                                    show_legend = FALSE)
                }
)


col_fun1 <- circlize::colorRamp2(c(0, 40), c("white", "#66C2A5"))
col_fun2 <- circlize::colorRamp2(c(0, 40), c("white", "#3E2789"))
col_fun3 <- circlize::colorRamp2(c(0, 30), c("white", "#CC075E"))
col_fun4 <- circlize::colorRamp2(c(0, 15), c("white", "dodgerblue"))
## 🔹 Draw heatmaps ----
hmaps_atac_junb_l <- map2(mtx_atac_junb_l, names(mtx_atac_junb_l),
                        ~ {
                          annot <- annot_junb_l[[.y]]
                          ComplexHeatmap::Heatmap(.x, show_heatmap_legend = FALSE,
                                                  col = col_fun2,
                                                  cluster_rows = FALSE, cluster_columns = FALSE, 
                                                  show_row_names = FALSE, show_column_names = FALSE,
                                                  row_order = order_for_heatmaps, border = TRUE,
                                                  top_annotation = annot,
                                                  column_title_gp = grid::gpar(fontsize = 7))
                        }
)


hmaps_jun_junb_wt_l <- map2(mtx_junb_junb_wt_l, names(mtx_junb_junb_wt_l),
                       ~ {
                         annot <- annot_junb_junb_l[[.y]]
                         ComplexHeatmap::Heatmap(.x, show_heatmap_legend = FALSE,
                                                 col = col_fun1,
                                                 cluster_rows = FALSE, cluster_columns = FALSE, 
                                                 show_row_names = FALSE, show_column_names = FALSE,
                                                 row_order = order_for_heatmaps, border = TRUE,
                                                 top_annotation = annot,
                                                 #column_title = titles_hmaps_wt[[name_mtx]], 
                                                 column_title_gp = grid::gpar(fontsize = 7))
                       }
)

hmaps_hac_junb_wt_l <- map2(mtx_h3ac_junb_wt_l, names(mtx_h3ac_junb_wt_l),
                       ~ {
                         annot <- annot_junb_junb_l[[.y]]
                         ComplexHeatmap::Heatmap(.x, show_heatmap_legend = FALSE,
                                                 col = col_fun3,
                                                 cluster_rows = FALSE, cluster_columns = FALSE, 
                                                 show_row_names = FALSE, show_column_names = FALSE,
                                                 row_order = order_for_heatmaps, border = TRUE,
                                                 top_annotation = annot,
                                                 #column_title = titles_hmaps_wt[[name_mtx]], 
                                                 column_title_gp = grid::gpar(fontsize = 7))
                       }
)

hmaps_hme_junb_wt_l <- map2(mtx_h3me_junb_wt_l, names(mtx_h3me_junb_wt_l),
                       ~ {
                         annot <- annot_junb_junb_l[[.y]]
                         ComplexHeatmap::Heatmap(.x, show_heatmap_legend = FALSE,
                                                 col = col_fun4,
                                                 cluster_rows = FALSE, cluster_columns = FALSE, 
                                                 show_row_names = FALSE, show_column_names = FALSE,
                                                 row_order = order_for_heatmaps, border = TRUE,
                                                 top_annotation = annot,
                                                 #column_title = titles_hmaps_wt[[name_mtx]], 
                                                 column_title_gp = grid::gpar(fontsize = 7))
                       }
)


partition_hp <- Heatmap(partition, col=structure(c("salmon2", "palevioletred"), names = c(paste("JuNB-independent\n", length(junb_indep)),
                                                                                          paste("JUNB-dependent\n", length(junb_dep)))), 
                        show_column_names = FALSE, name = " ", show_heatmap_legend = FALSE, 
                        row_order = order_for_heatmaps, row_title_gp = gpar(fontsize = 10), 
                        show_row_names = FALSE, width=unit(1,'mm'))

partition_hp2 <- Heatmap(partition, col=structure(c("salmon2", "palevioletred"), names = c(paste("JuNB-independent\n", length(junb_indep)),
                                                                                          paste("JUNB-dependent\n", length(junb_dep)))), 
                        show_column_names = FALSE, name = " ", show_heatmap_legend = FALSE, 
                        row_order = order_for_heatmaps, row_title = c("",""), row_title_gp = gpar(fontsize = 0), 
                        show_row_names = FALSE, width=unit(1,'mm'))


ht_list_atac_junb <- hmaps_atac_junb_l[["WT_V1"]] +
  hmaps_atac_junb_l[["KD_C1"]] 

ht_list_junb_junb <-   hmaps_jun_junb_wt_l[["junb_WT_VEH"]] +
  hmaps_jun_junb_wt_l[["junb_KD_CCM"]] 

ht_list_hk27_junb <-   hmaps_hac_junb_wt_l[["H3K27ac_WT_VEH"]] +
  hmaps_hac_junb_wt_l[["H3K27ac_KD_CCM"]] 

ht_list_hkme_junb <-   hmaps_hme_junb_wt_l[["H3K4me1_WT_VEH"]] +
  hmaps_hme_junb_wt_l[["H3K4me1_KD_CCM"]] 


## Save ATAC time points heatmap WT regions----
png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_junb_regions.png", width = 2.1, height = 3.5, units = "in", res = 300)
draw(partition_hp + ht_list_atac_junb,  #heatmap_legend_list = combined_legend_atac_wt, heatmap_legend_side = "bottom",
     gap = unit(1, "mm"), split = partition,
     column_title = "          Chromatin\n          accessibility", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/junb_junb_regions.png", width = 1.84, height = 3.5, units = "in", res = 300)
draw(partition_hp2 + ht_list_junb_junb,  #heatmap_legend_list = combined_legend_atac_wt, heatmap_legend_side = "bottom",
     gap = unit(1, "mm"), split = partition,
     column_title = "   JUNB\n    ChIP", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/k27a_junb_regions.png", width = 1.84, height = 3.5, units = "in", res = 300)
draw(partition_hp2 + ht_list_hk27_junb,  #heatmap_legend_list = combined_legend_atac_wt, heatmap_legend_side = "bottom",
     gap = unit(1, "mm"), split = partition,
     column_title = "   H3K27ac\n    ChIP", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/k4me_junb_regions.png", width = 1.84, height = 3.5, units = "in", res = 300)
draw(partition_hp2 + ht_list_hkme_junb,  #heatmap_legend_list = combined_legend_atac_wt, heatmap_legend_side = "bottom",
     gap = unit(1, "mm"), split = partition,
     column_title = "   H3K4me1\n    ChIP", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

### Junb_target genes ----

file_nojunb <- read.table("../results/ATAC/data/target_genes_of_junb_negative_atac_cluster.tsv")
file_junb <- read.table("../results/ATAC/data/target_genes_of_junb_positive_atac_cluster.tsv")
# 
genes_junb <- file_junb$V7
genes_nojunb <- file_nojunb$V7


lcpm.rlog <- read.csv(file = "../results/Human_cells/D3_RNA/data/log2cpm.csv", row.names = 1) 
colnames(lcpm.rlog)
lcpm.rlog <- lcpm.rlog[c(1:6, 13:18, 7:12, 19:24)] %>% rownames_to_column("gene")

### ATAC junb dep and ind - using WT and KD CCM ----

lcpm.rlog_genes_nojunb <- data.frame(gene = genes_nojunb) %>% left_join(lcpm.rlog[c(1:4, 17:19)], by = "gene")
lcpm.rlog_genes_junb <- data.frame(gene = genes_junb) %>% left_join(lcpm.rlog[c(1:4, 17:19)], by = "gene")

# Add a column to indicate source
nojunb_df <- lcpm.rlog_genes_nojunb %>% mutate(GeneSet = "JUNB-independent")
junb_df <- lcpm.rlog_genes_junb %>% mutate(GeneSet = "JUNB-dependent")

# Combine the two datasets while keeping duplicates
merged_unfilt_df <- bind_rows(nojunb_df, junb_df)


# Ensure the order is respected in the summarization
gene_counts <- merged_unfilt_df %>%
  arrange(match(GeneSet, c("JUNB-independent", "JUNB-dependent"))) %>%  # Ensure order
  group_by(gene) %>%
  summarise(`Gene regulation` = case_when(n() > 1 ~ "Both",  
                                          any(GeneSet == "JUNB-dependent") & !any(GeneSet == "JUNB-independent") ~ "JUNB-positive",
                                          TRUE ~ "JUNB-negative"), 
            .groups = "drop") 

merged_df  <- merged_unfilt_df %>% 
  left_join(gene_counts, by  ="gene") %>% drop_na()

my_colors <- broman::crayons(c("Asparagus","indigo", "Timberwolf"))

annotation_colors <- list("Gene regulation" = c("JUNB-positive" = my_colors[[1]], "JUNB-negative" = my_colors[[2]], "Both" = my_colors[[3]]))


annot <- data.frame(Condition = str_replace_all(colnames(merged_df[-c(1,8:9)]), c("\\d" = "", "WT_V" = "Control", "KD_C" = "siPDCD10\nCCM-like env")),
                    row.names = colnames(merged_df[-c(1,8:9)]))

annot$Condition <- factor(annot$Condition, levels = unique(annot$Condition))

#colors
cond_colors <- color_pal[c(1,6)]
names(cond_colors) <- unique(str_replace_all(colnames(merged_df[-c(1,8:9)]), c("\\d" = "", "WT_V" = "Control", "KD_C" = "siPDCD10\nCCM-like env")))

my_colour <-  list("Condition" = cond_colors)

left_annotation <- rowAnnotation("Gene regulation" = merged_df$`Gene regulation`,
                                 col = annotation_colors,
                                 show_annotation_name = FALSE,  show_legend = FALSE)

top_annot_gTF <- HeatmapAnnotation(Condition = annot$Condition,
                                   col = my_colour,
                                   annotation_height = unit(2, "mm"),
                                   show_annotation_name = FALSE, show_legend = FALSE) 



gene_label_motif <- c("SERPINE1", "ICAM1", "TGFB2", "VEGFC", "CCL2", "VEGFA")
#Get the Row Names
rownames_genes_neg <- merged_df[merged_df$GeneSet == "JUNB-independent", "gene"]
rownames_genes_pos <- merged_df[merged_df$GeneSet == "JUNB-dependent", "gene"]

#Match TF Motifs to Row Names
genes_positions_neg <- match(gene_label_motif, rownames_genes_neg)
genes_positions_pos <- match(gene_label_motif, rownames_genes_pos)

genes_positions_anno <- c(genes_positions_neg[!is.na(genes_positions_neg)], genes_positions_pos[!is.na(genes_positions_pos)]+length(rownames_genes_neg))

labels_anno <- c(gene_label_motif[!is.na(genes_positions_neg)], gene_label_motif[!is.na(genes_positions_pos)])

#Add Annotation to Heatmap
ha_genes <- rowAnnotation(gene_labels = anno_mark(at = genes_positions_anno, 
                                                  labels = labels_anno, 
                                                  labels_gp = gpar(fontsize = 9, col = "black"),
                                                  padding = unit(1, "mm"),
                                                  link_width = unit(2, "mm")))

partition_genes <- fct_inorder(c(rep("JuNB-independent", length(rownames_genes_neg)),
                             rep("JUNB-dependent", length(rownames_genes_pos))))


partition_genes_hp <- Heatmap(partition_genes, col=structure(c("salmon2", "palevioletred"), names = c("JuNB-independent", "JUNB-dependent")), 
                              show_column_names = FALSE, name = " ", show_heatmap_legend = FALSE, 
                              row_title_gp = gpar(fontsize = 10), 
                              show_row_names = FALSE, width=unit(1,'mm'))

genes_for_hmap <- t(scale(t(merged_df[, 2:7])))

control_annot <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = color_pal[[1]], 
                                                col = "black")), 
                  height = unit(2.5, "mm"),
                  show_legend = FALSE)

siPDCD10_annot <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = color_pal[[6]], 
                                                                  col = "black")), 
                                    height = unit(2.5, "mm"),
                                    show_legend = FALSE)

col_fun_junb_rna <- circlize::colorRamp2(seq(-3, 3, length = 100), plasma(100))


control_hmap <- ComplexHeatmap::Heatmap(genes_for_hmap[ ,1:3], show_heatmap_legend = FALSE,
                                          col = col_fun_junb_rna, 
                                          cluster_rows = FALSE, cluster_columns = FALSE,
                                          show_row_names = FALSE, show_column_names = FALSE,
                                           border = TRUE, left_annotation = left_annotation,
                                          top_annotation = control_annot)

siPDCD10_hmap <- ComplexHeatmap::Heatmap(genes_for_hmap[ ,4:6], show_heatmap_legend = FALSE,
                                        col = col_fun_junb_rna, 
                                        cluster_rows = FALSE, cluster_columns = FALSE, 
                                        show_row_names = FALSE, show_column_names = FALSE,
                                        border = TRUE, right_annotation = ha_genes,
                                        top_annotation = siPDCD10_annot)


png("../results/Human_cells/timepoint_reg/heatmaps_coverage/genes_junb_regions.png", width = 2.9, height = 3.5, units = "in", res = 300)
draw(partition_genes_hp + control_hmap + siPDCD10_hmap,  #heatmap_legend_list = combined_legend_atac_wt, heatmap_legend_side = "bottom",
     gap = unit(1, "mm"), split = partition_genes,
     column_title = "Gene     \nExpression     ", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()


atac_junb_legend <- Legend(col_fun = col_fun2, 
                           title = "RPKM", 
                           at = c(0,  20, 40),  # Define ticks on the legend
                           labels = c("0",  "20", "40"),direction  = "horizontal",  title_position = "topcenter",
                           title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                           labels_gp = gpar(fontsize = 8, col = "black"),
                           legend_height = unit(2.5, "cm"),
                           legend_width = unit(2, "cm"),
                           title_gap = unit(1, "mm"))
jun_junb_legend <- Legend(col_fun = col_fun1, 
                          title = "RPKM", direction = "horizontal",
                          at = c(0,  20, 40),  # Define ticks on the legend
                          labels = c("0",  "20", "40"), title_position = "topcenter",
                          title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                          labels_gp = gpar(fontsize = 8, col = "black"),
                          legend_height = unit(2, "cm"),
                          legend_width = unit(2, "cm"),
                          title_gap = unit(1, "mm"))
hac_junb_legend <- Legend(col_fun = col_fun3, 
                          title = "RPKM", direction = "horizontal",
                          at = c(0,  15, 30),  # Define ticks on the legend
                          labels = c("0",  "15", "30"), title_position = "topcenter",
                          title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                          labels_gp = gpar(fontsize = 8, col = "black"),
                          legend_height = unit(2, "cm"),
                          legend_width = unit(2, "cm"),
                          title_gap = unit(1, "mm"))
hme_junb_legend <- Legend(col_fun = col_fun4, 
                          title = "RPKM", direction = "horizontal", title_position = "topcenter",
                          at = c(0,  8, 15),  # Define ticks on the legend
                          labels = c("0",  "8", "15"),
                          title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                          labels_gp = gpar(fontsize = 8, col = "black"),
                          legend_height = unit(2, "cm"),
                          legend_width = unit(2, "cm"),
                          title_gap = unit(1, "mm"))

cts_rna_junb_legend <- Legend(col_fun = col_fun_junb_rna, 
                              direction = "horizontal",
                         title = "Z-score log2 (CPM +1)", title_position = "topcenter",
                         title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                         labels_gp = gpar(fontsize = 8, col = "black"),
                         legend_height = unit(2, "cm"),
                         legend_width = unit(2, "cm"),
                         title_gap = unit(1, "mm"))


cond_legend <- Legend(labels = c("Control\n", "siPDCD10\nCCM like env"),
                      title = "Condition", 
                      row_gap = unit(2, "mm"),
                      legend_gp = gpar(fill = color_pal[c(1, 6)]),
                      title_gp = gpar(fontsize = 11, fontface = "bold"),   # Legend title appearance
                      labels_gp = gpar(fontsize = 10, col = "black"),
                      title_gap = unit(1, "mm"))
junbreg_legend <- Legend(labels = c("JUNB-positive", "JUNB-negative", "Both"),
                      title = "Gene regulation", 
                      row_gap = unit(2, "mm"),
                      legend_gp = gpar(fill = my_colors),
                      title_gp = gpar(fontsize = 11, fontface = "bold"),   # Legend title appearance
                      labels_gp = gpar(fontsize = 10, col = "black"),
                      title_gap = unit(1, "mm"))


annotation_colors <- list("Gene regulation" = c("JUNB-positive" = "violet", "JUNB-negative" = "blue", "Both" = "gray"))


png(filename = "../results/Human_cells/timepoint_reg/heatmaps_coverage/legend_hu_junbregions.png", width = 1.2, height = 3.5, units = "in", res = 300)
ComplexHeatmap::draw(packLegend(cond_legend, junbreg_legend, direction = "vertical"))
dev.off()


combined_legend_junbregions <- packLegend(atac_junb_legend, jun_junb_legend, hac_junb_legend, hme_junb_legend, cts_rna_junb_legend, 
                                 direction = "horizontal", gap = unit(0.45, "in"))

png(filename = "../results/Human_cells/timepoint_reg/heatmaps_coverage/legend_junbregions.png", width = 6.8, height = 0.5, units = "in", res = 300)
ComplexHeatmap::draw(combined_legend_junbregions, x = unit(0.5, "npc"))
dev.off()


# Optionally save session info
sink("../logs/sessioninfo_human_15_coverage_junb_dep_indep_regions.txt")
sessionInfo()
sink()