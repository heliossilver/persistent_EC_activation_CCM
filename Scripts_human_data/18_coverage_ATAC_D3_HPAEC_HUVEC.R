
pacman::p_load("dplyr", "tibble", "tidyr", "ggplot2", "readr",
               "rtracklayer", "magick", "GenomicRanges", "purrr",
               "vroom", "stringr", "viridis", "RColorBrewer", 
               "glue", "EnrichedHeatmap", "forcats", "enrichR")

dir.create("../results/Human_cells/timepoint_reg/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("../results/Human_cells/timepoint_reg/heatmaps_coverage", showWarnings = FALSE, recursive = TRUE)

color_pal <- RColorBrewer::brewer.pal(12, "Paired")

# color_pal_map <- setNames(c(rep(color_pal[1:2], each = 3), rep(color_pal[1:2], each = 3), 
#                             rep(color_pal[5:6], each = 3), rep(color_pal[5:6], each = 3)),
#                           condition_names)

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

##human data HPAEC, HUVEC, hCMEC----
hpaec_ccm <- read_tsv("../results/Human_cells/common/venn/A_unique_ccm.bed", col_names = c("seqnames", "start", "end")) %>% 
  makeGRangesFromDataFrame()
huvec_ccm <- read_tsv("../results/Human_cells/common/venn/B_unique_ccm.bed", col_names = c("seqnames", "start", "end")) %>% 
  makeGRangesFromDataFrame()
hcmec_ccm <- read_tsv("../results/Human_cells/common/venn/C_unique_ccm.bed", col_names = c("seqnames", "start", "end")) %>% 
  makeGRangesFromDataFrame()
common_ec_ccm <- read_tsv("../results/Human_cells/common/venn/triple_ccm_fromA.bed", col_names = c("seqnames", "start", "end")) %>% 
  makeGRangesFromDataFrame()

hpaec_veh <- read_tsv("../results/Human_cells/common/venn/A_unique_vehicle.bed", col_names = c("seqnames", "start", "end")) %>% 
  makeGRangesFromDataFrame()
huvec_veh <- read_tsv("../results/Human_cells/common/venn/B_unique_vehicle.bed", col_names = c("seqnames", "start", "end")) %>% 
  makeGRangesFromDataFrame()
hcmec_veh <- read_tsv("../results/Human_cells/common/venn/C_unique_vehicle.bed", col_names = c("seqnames", "start", "end")) %>% 
  makeGRangesFromDataFrame()
common_ec_veh <- read_tsv("../results/Human_cells/common/venn/triple_vehicle_fromA.bed", col_names = c("seqnames", "start", "end")) %>% 
  makeGRangesFromDataFrame()


A <- read.csv("../results/Human_cells/HPAEC_ATAC/edgeR/deg_difexp_tables/HPAEC_ccm_vs_HPAEC_veh_diffacc.csv") %>% 
  dplyr:::filter(diffaccessible == "Up") %>% 
  separate(peaks, into = c("seqnames", "start", "end"), sep = "-") %>% makeGRangesFromDataFrame()
B <- read.csv("../results/Human_cells/HuvEC_ATAC/edgeR/deg_difexp_tables/HUVEC_ccm_vs_HUVEC_veh_diffacc.csv") %>% 
  dplyr:::filter(diffaccessible == "Up") %>% 
  separate(peaks, into = c("seqnames", "start", "end"), sep = "-")%>% makeGRangesFromDataFrame()
C <- read_table("../results/Human_cells/data/WT_CCM_vs_WT_VEH_peaks.txt", 
                col_names = c("seqnames", "start", "end")) %>% makeGRangesFromDataFrame()

ABC <- c(A, B, C)

hpaec_veh <- hpaec_veh[!hpaec_veh %over% A]
huvec_veh <- huvec_veh[!huvec_veh %over% B]
hcmec_veh <- hcmec_veh[!hcmec_veh %over% C]
common_ec_veh <- common_ec_veh[!common_ec_veh %over% ABC]


# --- 1) Put your four GRanges into a named list -----

hpaec_ccm_strict <- hpaec_ccm[ !overlapsAny(hpaec_ccm, hpaec_veh) ]
huvec_ccm_strict <- huvec_ccm[ !overlapsAny(huvec_ccm, huvec_veh) ]
hcmec_ccm_strict <- hcmec_ccm[ !overlapsAny(hcmec_ccm, hcmec_veh) ]
common_ccm_strict <- common_ec_ccm[ !overlapsAny(common_ec_ccm, common_ec_veh) ]

gr_sets_ccm_human <- c(hpaec_ccm_strict,  
                       huvec_ccm_strict,
                       hcmec_ccm_strict,
                       common_ccm_strict)  

partition_hu_treated <- fct_inorder(c(rep(paste("HPAEC\n", "n =", length(hpaec_ccm_strict)), length(hpaec_ccm_strict)),
                                      rep(paste("HUVEC\n", "n =", length(huvec_ccm_strict)), length(huvec_ccm_strict)),
                                      rep(paste("hCMEC/D3\n", "n =", length(hcmec_ccm_strict)),length(hcmec_ccm_strict)),
                                      rep(paste("Common\n", "n =", length(common_ccm_strict)), length(common_ccm_strict))))
length(partition_hu_treated)

gene_enh_pairs <- read_csv("../results/Human_cells/D3_ATAC/data/pcc_ccre_genes_500kb_filt.csv")

gene_enh_pairs_gr <- gene_enh_pairs %>% 
  separate(peaks, into = c("chr", "start", "end"), sep = "-", remove = FALSE) %>% 
  dplyr::select(chr, start, end, peaks, gene) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

gr_sets_ccm_human_w_metadata <- gr_sets_ccm_human
mcols(gr_sets_ccm_human_w_metadata)$classification <- c(rep("HPAEC", length(hpaec_ccm_strict)),
                                             rep("HUVEC",  length(huvec_ccm_strict)),
                                             rep("hCMEC/D3", length(hcmec_ccm_strict)),
                                             rep("Common",         length(common_ccm_strict)))
peaks_ids <- paste(seqnames(gr_sets_ccm_human), 
                   start(gr_sets_ccm_human), 
                   end(gr_sets_ccm_human), sep = "-")

# 1. Find overlaps between the two GRanges objects
hits <- findOverlaps(gr_sets_ccm_human_w_metadata,
                     gene_enh_pairs_gr)

# 2. For each region in gr_sets_ccm_human_w_metadata that hit,
#    collect the gene(s) it's associated with.
qry_idx <- queryHits(hits)
sub_idx <- subjectHits(hits)

# pull the gene column (rename if yours is different)
gene_vec <- mcols(gene_enh_pairs_gr)$gene[sub_idx]
# tapply(gene_vec, qry_idx, function(x) paste(unique(x), collapse = ";"))
# 3. Keep only the ranges that had a match
ccm_matched <- gr_sets_ccm_human_w_metadata[(qry_idx)]

# 4. Add the gene annotation as a new metadata column
#mcols(ccm_matched)$linked_genes <- genes_by_query[ match(seq_along(ccm_matched), unique(qry_idx)) ]
mcols(ccm_matched)$linked_gene <- mcols(gene_enh_pairs_gr)$gene[sub_idx]
write.csv(as.data.frame(ccm_matched), "../results/Human_cells/timepoint_reg/tables/linked_genes_commonCCM_humanEC.csv", row.names = FALSE)


##bw 
bw_file <- list.files("../results/Human_cells/D3_ATAC/bw/", full.names = T)
key_bw_names <- read.csv("../results/Human_cells/D3_ATAC/data/key_atac.csv", header = FALSE) 


bw_true_names <- setNames(key_bw_names$V1, key_bw_names$V2)

bw_names <- str_replace_all(basename(bw_file), c(bw_true_names, "_sorted_rmdup.bw" = ""))

bw_list_hu_comp <- setNames(lapply(bw_file[c(37:44, 29, 30, 25, 33, 34, 27)], import), 
                            bw_names[c(37:44, 29, 30, 25, 33, 34, 27)])

norm_atac_mtx_hu_l <- normalize_bw_list(bw_list_hu_comp, gr_sets_ccm_human)
#saveRDS(norm_atac_mtx_hu_l, "../results/Human_cells/bw_norm_mtx/norm_atac_mtx_hu_l.rds")
#norm_atac_mtx_hu_l <- readRDS("../results/Human_cells/bw_norm_mtx/norm_atac_mtx_hu_l.rds")

norm_atac_mtx_hu_l

# Compute quantiles for norm----
quantiles_atac_hu <- lapply(norm_atac_mtx_hu_l, compute_quantiles)
quantiles_atac_hu

col_fun2 <- circlize::colorRamp2(c(0, 50), c("white", "#3E2789"))

conditions_hu <- c(rep(c("HPAEC\nVehicle", "HPAEC\nTNF/DMOG"), each =2),
                   rep(c("HUVEC\nVehicle", "HUVEC\nTNF/DMOG"), each =2),
                   rep(c("hCMEC\nVehicle", "hCMEC\nTNF/DMOG"), each =2))

titles_hmaps_hu <- setNames(conditions_hu,
                            names(norm_atac_mtx_hu_l))


mtx_atac_hu_l <- map2(norm_atac_mtx_hu_l, names(norm_atac_mtx_hu_l), 
                      ~{
                        mtx <- .x %>%
                          as.data.frame() %>% 
                          as.matrix()
                        rownames(mtx) <- peaks_ids
                        mtx
                      })

# 2) One score per row by averaging rowMeans across all matrices
row_scores <- map_dbl(peaks_ids, function(id) {
  vals <- map_dbl(mtx_atac_hu_l, ~ mean(as.numeric(.x[id, ]), na.rm = TRUE))
  mean(vals, na.rm = TRUE)  # global score for this row
})

# 3) Build a single ordering (within class, desc by score)
ord_ids <- tibble(id = peaks_ids,
                  class = partition_hu_treated,
                  score = row_scores) %>%
  group_by(class) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  pull(id)

# 4) Reorder all matrices using the SAME order
mtx_atac_hu_l <- map2(mtx_atac_hu_l, names(mtx_atac_hu_l), ~{
  .x[ord_ids, , drop = FALSE]
})


color_pal <- RColorBrewer::brewer.pal(12, "Paired")

# color_pal_map <- setNames(c(rep(c(color_pal[11],"yellow3") , each = 2), 
#                             rep(color_pal[7:8], each = 2),
#                             rep(color_pal[1:2], each = 2)),
#                           conditions_hu)

color_pal_map <- setNames(rep(rep(c("#B2DF8A", "#33A02C"), each = 2), 3),
                          conditions_hu)



annot_l <- map2(titles_hmaps_hu, names(titles_hmaps_hu), 
                ~{
                  HeatmapAnnotation(Source = anno_block(gp = gpar(fill = color_pal_map[[.x]], 
                                                                  col = "black")), 
                                    height = unit(2.5, "mm"),
                                    show_legend = FALSE)
                }
)

mtx_example <- mtx_atac_hu_l[["WT_C1"]]
genes_peaks_associted <- read.csv("../results/Human_cells/timepoint_reg/tables/linked_genes_commonCCM_humanEC.csv")  %>% 
  unite(peaks, c(seqnames, start, end), sep = "-")

length(rownames(mtx_example))
# Peaks linked to SERPINE1
peak_gene_ids <- genes_peaks_associted %>%
  dplyr::filter(classification == "Common",
                linked_gene %in% c("SERPINE1", "VCAM1", "ICAM1", "CCL2")) %>% 
  mutate(position =  match(peaks, rownames(mtx_example))) %>% 
  arrange(desc(position)) %>% distinct(linked_gene, .keep_all = T)



ha <- rowAnnotation(gene_label = anno_mark(at = peak_gene_ids$position,
                                           labels = peak_gene_ids$linked_gene,
                                           labels_gp = gpar(fontsize = 5, col = "black"),
                                           padding = unit(1, "mm"),
                                           link_width = unit(2, "mm")))


### drawing heatmaps----
hmaps_atac_hu_l <- map2(mtx_atac_hu_l, names(mtx_atac_hu_l),
                        ~ {
                          annot <- annot_l[[.y]]
                          ComplexHeatmap::Heatmap(.x, show_heatmap_legend = FALSE,
                                                  col = col_fun2,
                                                  cluster_rows = FALSE, cluster_columns = FALSE, 
                                                  show_row_names = FALSE, show_column_names = FALSE,
                                                  #row_order = order_for_heatmaps_hu, 
                                                  border = TRUE,
                                                  top_annotation = annot,
                                                  #column_title = titles_hmaps_wt[[name_mtx]], 
                                                  column_title_gp = grid::gpar(fontsize = 7))
                        }
)

wtc1_hmap <- ComplexHeatmap::Heatmap(mtx_atac_hu_l[["WT_C1"]], show_heatmap_legend = FALSE,
                        col = col_fun2,
                        cluster_rows = FALSE, cluster_columns = FALSE, 
                        show_row_names = FALSE, show_column_names = FALSE,
                        #row_order = order_for_heatmaps_hu, 
                        border = TRUE,
                        top_annotation = annot_l[["WT_C1"]],
                        right_annotation = ha,
                        #column_title = titles_hmaps_wt[[name_mtx]], 
                        column_title_gp = grid::gpar(fontsize = 7))

partition_hp_hu <- Heatmap(partition_hu_treated, col=structure(c("#BAD1CD", "#F2D1C9", "#E086D3",  "#64113F"), 
                                                               names = c(paste("HPAEC\n", "n =", length(hpaec_ccm_strict)), 
                                                                         paste("HUVEC\n", "n =", length(huvec_ccm_strict)),  
                                                                         paste("hCMEC/D3\n", "n =", length(hcmec_ccm_strict)),
                                                                         paste("Common\n", "n =", length(common_ccm_strict)))), 
                           show_column_names = FALSE, name = " ", show_heatmap_legend = FALSE, 
                           #row_order = order_for_heatmaps_hu,   
                           row_title_gp = gpar(fontsize = 7),
                           show_row_names = FALSE, width=unit(1,'mm'))

partition_hp_hu2 <- Heatmap(partition_hu_treated, col=structure(c("#BAD1CD", "#F2D1C9", "#E086D3",  "#64113F"), 
                                                               names = c(paste("HPAEC\n", "n =", length(hpaec_ccm_strict)), 
                                                                         paste("HUVEC\n", "n =", length(huvec_ccm_strict)),  
                                                                         paste("hCMEC/D3\n","n =",  length(hcmec_ccm_strict)),
                                                                         paste("Common\n","n =",  length(common_ccm_strict)))), 
                           show_column_names = FALSE, name = " ", show_heatmap_legend = FALSE, 
                           #row_order = order_for_heatmaps_hu,  
                           row_title = c("","","", ""),
                           row_title_gp = gpar(fontsize = 0),
                           show_row_names = FALSE, width=unit(1,'mm'))

atac_legend <- Legend(col_fun = col_fun2, 
                      title = "RPKM", title_position = "topcenter",
                      at = c(0,  25, 50),  # Define ticks on the legend
                      labels = c("0",  "25", "50"), direction = "horizontal",
                      title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                      labels_gp = gpar(fontsize = 8, col = "black"),
                      legend_height = unit(1, "cm"),
                      legend_width = unit(2.5, "cm"), 
                      title_gap = unit(1, "mm"))


combined_legend <- packLegend(atac_legend, direction = "horizontal")

ht_list_hpaec <- hmaps_atac_hu_l[["HPAEC_V1"]] +
  hmaps_atac_hu_l[["HPAEC_C1"]] 
ht_list_huvec <-   hmaps_atac_hu_l[["HUVEC_V1"]] +
  hmaps_atac_hu_l[["HUVEC_C1"]] 
ht_list_hcmec <- hmaps_atac_hu_l[["WT_V1"]] +
  wtc1_hmap


png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_hpaec.png", width = 1.6, height = 3, units = "in", res = 300)
draw(partition_hp_hu + ht_list_hpaec, #heatmap_legend_list = combined_legend, 
     gap = unit(1, "mm"), split = partition_hu_treated)
dev.off()
png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_huvec.png", width = 1.43, height = 3, units = "in", res = 300)
draw(partition_hp_hu2 + ht_list_huvec, #heatmap_legend_list = combined_legend, .1
     gap = unit(1, "mm"), split = partition_hu_treated)
dev.off()
png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_hcmec.png", width = 1.9, height = 3, units = "in", res = 300)
draw(partition_hp_hu2 + ht_list_hcmec, 
     gap = unit(1, "mm"), split = partition_hu_treated)
dev.off()


cond_legend_hu <- Legend(labels = c("Vehicle\n", "CCM-like\nenvironment"),
                         title = "Condition", 
                         row_gap = unit(2, "mm"),
                         legend_gp = gpar(fill = unique(color_pal_map)[1:2]),
                         title_gp = gpar(fontsize = 9, fontface = "bold"),   # Legend title appearance
                         labels_gp = gpar(fontsize = 8.5, col = "black"),
                         title_gap = unit(1, "mm"))

png(filename = "../results/Human_cells/timepoint_reg/heatmaps_coverage/legend_hu.png", width = 1.43, height = 0.5, units = "in", res = 300)
ComplexHeatmap::draw(combined_legend)
dev.off()
png(filename = "../results/Human_cells/timepoint_reg/heatmaps_coverage/legend_hu_conditions.png", width = 1, height = 3.5, units = "in", res = 300)
ComplexHeatmap::draw(packLegend(cond_legend_hu, direction = "vertical"))
dev.off()


# Optionally save session info
sink("../logs/sessioninfo_human_18_coverage_ATAC_D3_HPAEC_HUVEC.txt")
sessionInfo()
sink()

