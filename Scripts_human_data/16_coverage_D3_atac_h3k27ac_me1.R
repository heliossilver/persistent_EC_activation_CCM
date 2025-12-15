
pacman::p_load("dplyr", "tibble", "tidyr", 
               "ggplot2", "readr",
               "rtracklayer", "magick", 
               "GenomicRanges", "purrr",
               "vroom", "stringr", "viridis",
               "RColorBrewer", "glue", 
               "EnrichedHeatmap", "forcats", "enrichR")

dir.create("../results/Human_cells/timepoint_reg/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("../results/Human_cells/timepoint_reg/heatmaps_coverage", showWarnings = FALSE, recursive = TRUE)

###Functions ----
# Helper to read an edgeR/DE table -> GRanges of UP peaks--
read_up_as_gr <- function(csv_path, peak_col = "peaks", direction = c("Up", "Down")) {
  df <- read.csv(csv_path, check.names = FALSE) 
  
  logfc_col <- grep("^logFC", names(df), value = TRUE)
  fdr_col <- grep("^FDR", names(df), value = TRUE)
  
  df <- df %>%
    dplyr::rename(logFC = all_of(logfc_col), FDR = all_of(fdr_col)) %>%
    {
      if (direction == "Up") {
        dplyr::filter(., logFC >= 0.5,
                      FDR < 0.05)
      } else {
        dplyr::filter(., logFC <= -0.5,
                      FDR < 0.05) 
      }  
    } %>%                     
    tidyr::separate({{peak_col}}, c("chr","start","end"), sep = "-", remove = FALSE) %>%
    mutate(start = as.integer(start),
           end   = as.integer(end)) %>%
    distinct(chr, start, end, .keep_all = TRUE)
  
  makeGRangesFromDataFrame(df, keep.extra.columns = TRUE,
                           seqnames.field = "chr", start.field = "start", end.field = "end")
}

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
# saving bed files for motif enrichement--
gr_to_bed0 <- function(gr, file, extra_cols = NULL) {
  df <- as.data.frame(gr)
  
  # BED must be: chr, start(0-based), end
  bed <- data.frame(seqnames = df$seqnames,
                    start    = df$start - 1L,   # convert to 0-based
                    end      = df$end)
  
  # if you want to keep additional metadata (like score, name, strand, etc.)
  if (!is.null(extra_cols)) {
    bed <- cbind(bed, df[ , extra_cols, drop = FALSE])
  }
  
  # write to file
  write_tsv(bed, file, col_names = FALSE, quote = "none")
}

compute_quantiles <- function(mtx) {
  quantile(mtx, c(0.97, 0.99), na.rm = TRUE)
}


# ⚠️ Use the *within-KD* contrasts if your question is “opens at 1h/36h vs vehicle”----
# e.g., KD_T1h_vs_KD_Veh_diffexp.csv and KD_T36h_vs_KD_Veh_diffexp.csv
KD_1h_up  <- read_up_as_gr("../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables/K1_CCM_vs_W1_VEH_diffacc_full.csv", direction = "Up")
KD_36h_VEH_up <- read_up_as_gr("../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables/KD_CCM_vs_WT_VEH_diffacc_full.csv", direction = "Up")
KD_36h_1h_down <- read_up_as_gr("../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables/KD_CCM_vs_K1_CCM_diffacc_full.csv", direction = "Down")
KD_36h_1h_up <- read_up_as_gr("../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables/KD_CCM_vs_K1_CCM_diffacc_full.csv", direction = "Up")

WT_1h_up  <- read_up_as_gr("../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables/W1_CCM_vs_W1_VEH_diffacc_full.csv", direction = "Up")
WT_36h_VEH_up <- read_up_as_gr("../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables/WT_CCM_vs_WT_Veh_diffacc_full.csv", direction = "Up")
WT_36h_1h_down <- read_up_as_gr("../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables/WT_CCM_vs_W1_CCM_diffacc_full.csv", direction = "Down")
WT_36h_1h_up <- read_up_as_gr("../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables/WT_CCM_vs_W1_CCM_diffacc_full.csv", direction = "Up")


# ---- Primary membership (keeps whole peaks) WT ----
persistent_WT <- WT_1h_up[  overlapsAny(WT_1h_up,  WT_36h_VEH_up) ]
early_only_WT <- WT_1h_up[ !overlapsAny(WT_1h_up,  WT_36h_VEH_up) ]
late_only_WT  <- WT_36h_VEH_up[ !overlapsAny(WT_36h_VEH_up, WT_1h_up) ]

# ---- Strict confirmations with 36h vs 1h contrasts ----
# Early-only should *decrease* from 1h to 36h (36h vs 1h = Down)
early_only_WT_strict <- early_only_WT[ overlapsAny(early_only_WT, WT_36h_1h_down) ]

# Persistent should *not* be down from 1h to 36h; optionally also require it's Up (or at least not Down)
persistent_WT_strict <- persistent_WT[!overlapsAny(persistent_WT, WT_36h_1h_down) | overlapsAny(persistent_WT, WT_36h_1h_up)]

# Late-only should *increase* from 1h to 36h (36h vs 1h = Up)

late_only_WT_strict  <- late_only_WT[  overlapsAny(late_only_WT,  WT_36h_1h_up) ]

# --- 1) Put your four GRanges into a named list ----
gr_sets_wt <- c(early_only_WT_strict,  
                persistent_WT_strict,
                late_only_WT_strict)  

partition_names <- c(paste("Early response\nn =", c(length(early_only_WT_strict) + length(persistent_WT_strict))),
                     paste("Late response\n n =", length(late_only_WT_strict)))
                         
                         
partition_wt <- fct_inorder(rep(partition_names, c(c(length(early_only_WT_strict) + length(persistent_WT_strict), length(late_only_WT_strict)))))
length(partition_wt)


##bw wt----
bw_file <- list.files("../results/Human_cells/D3_ATAC/bw/", full.names = T)
key_bw_names <- read.csv("../results/Human_cells/D3_ATAC/data/key_atac.csv", header = FALSE) 


bw_true_names <- setNames(key_bw_names$V1, key_bw_names$V2)

bw_names <- str_replace_all(basename(bw_file), c(bw_true_names, "_sorted_rmdup.bw" = ""))

# names(bw_file) <- bw_names

bw_list <- setNames(lapply(bw_file[c(1:3, 4:6, 29:30, 25, 33:34, 27, 7:9, 10:12, 31:32, 26, 35:36, 28)], import), 
                    bw_names[c(1:3, 4:6, 29:30, 25, 33:34, 27, 7:9, 10:12, 31:32, 26, 35:36, 28)])


# 🔹 Normalize and save with correct filenames ----
norm_atac_mtx_wt_l <- normalize_bw_list(bw_list, gr_sets_wt) 

#saveRDS(norm_atac_mtx_wt_l, "../results/Human_cells/bw_norm_mtx/norm_atac_mtx_wt_l.rds")
#norm_atac_mtx_wt_l <- readRDS("../results/Human_cells/bw_norm_mtx/norm_atac_mtx_wt_l.rds")

# Compute quantiles for norm----
quantiles_atac_wt <- lapply(norm_atac_mtx_wt_l, compute_quantiles)
quantiles_atac_wt

col_fun2 <- circlize::colorRamp2(c(0, 90), c("white", "#3E2789"))
#col_fun2 <- circlize::colorRamp2(seq(0, 90, length = 100), magma(n = 90))


condition_names <- rep(c("Control\n(1h)", "Control\nCCM-like env (1h)", 
                         "Control\n(36h)", "Control\nCCM-like env (36h)",
                         "siPDCD10\n(1h)", "siPDCD10\nCCM-like env (1h)", 
                         "siPDCD10\n(36h)", "siPDCD10\nCCM-like env (36h)"), each = 3)

titles_hmaps_wt <- setNames(condition_names,
                            names(norm_atac_mtx_wt_l))


mtx_atac_wt_l <- map2(norm_atac_mtx_wt_l, names(norm_atac_mtx_wt_l), 
                      ~{
                        .x %>%
                          as.data.frame() %>%
                          as.matrix()
                      }
)


ht <- EnrichedHeatmap(norm_atac_mtx_wt_l[["WT_C1"]], pos_line = FALSE, show_heatmap_legend = FALSE,
                      col = col_fun2, axis_name = "", split = partition_wt,
                      top_annotation = NULL)
row_order_ht <- row_order(ht)

order_for_heatmaps <- as.vector(unlist(row_order_ht))

color_pal <- RColorBrewer::brewer.pal(12, "Paired")

color_pal_map <- setNames(c(rep(color_pal[1:2], each = 3), rep(color_pal[1:2], each = 3), 
                            rep(color_pal[5:6], each = 3), rep(color_pal[5:6], each = 3)),
                          condition_names)


annot_l <- map2(titles_hmaps_wt, names(titles_hmaps_wt), 
                ~{
                  HeatmapAnnotation(Source = anno_block(gp = gpar(fill = color_pal_map[[.x]], 
                                                                  col = "black")), 
                                    height = unit(2.5, "mm"),
                                    show_legend = FALSE)
                }
)


## 🔹 Draw heatmaps ----
hmaps_atac_wt_l <- map2(mtx_atac_wt_l, names(mtx_atac_wt_l),
                        ~ {
                          annot <- annot_l[[.y]]
                          ComplexHeatmap::Heatmap(.x, show_heatmap_legend = FALSE,
                                                  col = col_fun2,
                                                  cluster_rows = FALSE, cluster_columns = FALSE, 
                                                  show_row_names = FALSE, show_column_names = FALSE,
                                                  row_order = order_for_heatmaps, border = TRUE,
                                                  top_annotation = annot,
                                                  #column_title = titles_hmaps_wt[[name_mtx]], 
                                                  column_title_gp = grid::gpar(fontsize = 7))
                        }
)

partition_hp_wt <- Heatmap(partition_wt, col=structure(c("#EABFCB", 
                                                         #"#C191A1", 
                                                         "#A4508B"), names = partition_names), 
                           show_column_names = FALSE, name = " ", 
                           show_heatmap_legend = FALSE,  
                           row_order = order_for_heatmaps, row_title_gp = gpar(fontsize = 7), 
                           show_row_names = FALSE, width = unit(1,'mm'))

partition_hp_wt2 <- Heatmap(partition_wt, col=structure(c("#EABFCB", 
                                                          #"#C191A1", 
                                                          "#A4508B"), names = partition_names), 
                            show_column_names = FALSE, name = " ", 
                            show_heatmap_legend = FALSE,  row_title = c("",""),
                            row_order = order_for_heatmaps, row_title_gp = gpar(fontsize = 0), 
                            show_row_names = FALSE, width=unit(1,'mm'))

wt_legend <- Legend(col_fun = col_fun2, 
                    title = "RPKM", 
                    at = c(0,  45, 90),  # Define ticks on the legend
                    labels = c("0",  "45", "90"),direction  = "horizontal",  title_position = "topcenter",
                    title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                    labels_gp = gpar(fontsize = 8, col = "black"),
                    legend_height = unit(2.5, "cm"),
                    title_gap = unit(1, "mm"))

cond_legend_wt <- Legend(labels = c("Control\n", "Control\nCCM-like env"),
                         title = "Condition", 
                         row_gap = unit(2, "mm"),
                         legend_gp = gpar(fill = unique(color_pal_map)[1:2]),
                         title_gp = gpar(fontsize = 9, fontface = "bold"),   # Legend title appearance
                         labels_gp = gpar(fontsize = 8.5, col = "black"),
                         title_gap = unit(1, "mm"))
cond_legend_kd <- Legend(labels = c("siPDCD10\n", "siPDCD10\nCCM-like env"), 
                         title = "Condition",  
                         row_gap = unit(2, "mm"),
                         legend_gp = gpar(fill = unique(color_pal_map)[3:4]),
                         title_gp = gpar(fontsize = 9, fontface = "bold"),   # Legend title appearance
                         labels_gp = gpar(fontsize = 8.5, col = "black"),
                         title_gap = unit(1, "mm"))

combined_legend_wt <- packLegend(cond_legend_wt, direction = "vertical")
combined_legend_kd <- packLegend(cond_legend_kd, direction = "vertical")



ht_list_wt_wtregions_1h <- hmaps_atac_wt_l[["W1_V1"]] +
  hmaps_atac_wt_l[["W1_C1"]] 
ht_list_wt_wtregions_36h <-   hmaps_atac_wt_l[["WT_V1"]] +
  hmaps_atac_wt_l[["WT_C1"]] 

ht_list_kd_wtregions_1h <- hmaps_atac_wt_l[["K1_V1"]] +
  hmaps_atac_wt_l[["K1_C1"]] 

ht_list_kd_wtregions_36h <-   hmaps_atac_wt_l[["KD_V1"]] +
  hmaps_atac_wt_l[["KD_C1"]] 


## Save ATAC time points heatmap WT regions----
png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_timepoints_wt_wtregions1h.png", width = 1.45, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + ht_list_wt_wtregions_1h,  #heatmap_legend_list = wt_legend, heatmap_legend_side = "bottom",
     gap = unit(1, "mm"), split = partition_wt,
     column_title = "         1 hour", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()
png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_timepoints_wt_wtregions36h.png", width = 1.25, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt2 + ht_list_wt_wtregions_36h, #heatmap_legend_list = combined_legend_atac_wt, 
     show_heatmap_legend = FALSE, gap = unit(1, "mm"), split = partition_wt,
     column_title = "     36 hours", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_timepoints_kd_wtregions1h.png", width = 1.45, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + ht_list_kd_wtregions_1h,  #heatmap_legend_list = combined_legend_atac_wt, 
     show_heatmap_legend = FALSE, gap = unit(1, "mm"), split = partition_wt,
     column_title = "         1 hour", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()
png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_timepoints_kd_wtregions36h.png", width = 1.25, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt2 + ht_list_kd_wtregions_36h,# heatmap_legend_list = combined_legend_atac_wt, 
     show_heatmap_legend = FALSE, gap = unit(1, "mm"), split = partition_wt,
     column_title = "     36 hours", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png(filename = "../results/Human_cells/timepoint_reg/heatmaps_coverage/legendwt.png", width = 1.1, height = 3.14, units = "in", res = 300)
ComplexHeatmap::draw(combined_legend_wt)
dev.off()
png(filename = "../results/Human_cells/timepoint_reg/heatmaps_coverage/legendkd.png", width = 1.1, height = 3.14, units = "in", res = 300)
ComplexHeatmap::draw(combined_legend_kd)
dev.off()

## Save ATAC time points heatmap replicates----

titles_hmaps_rep <- rep(c("Rep. #1", "Rep. #2", "Rep. #3"), c(8))
names(titles_hmaps_rep) <- names(titles_hmaps_wt)

hmaps_atac_wt_rep_l <- map2(mtx_atac_wt_l, names(mtx_atac_wt_l),
                        ~ {
                          annot <- annot_l[[.y]]
                          ComplexHeatmap::Heatmap(.x, show_heatmap_legend = FALSE,
                                                  col = col_fun2,
                                                  cluster_rows = FALSE, cluster_columns = FALSE, 
                                                  show_row_names = FALSE, show_column_names = FALSE,
                                                  row_order = order_for_heatmaps, border = TRUE,
                                                  top_annotation = annot,
                                                  column_title = titles_hmaps_rep[[.y]], 
                                                  column_title_gp = grid::gpar(fontsize = 8))
                        }
)

ht_list_W1_V <- hmaps_atac_wt_rep_l[["W1_V1"]] + hmaps_atac_wt_rep_l[["W1_V2"]] + hmaps_atac_wt_rep_l[["W1_V3"]]
ht_list_WT_V <- hmaps_atac_wt_rep_l[["WT_V1"]] + hmaps_atac_wt_rep_l[["WT_V2"]] + hmaps_atac_wt_rep_l[["WT_V3"]]
ht_list_W1_C <- hmaps_atac_wt_rep_l[["W1_C1"]] + hmaps_atac_wt_rep_l[["W1_C2"]] + hmaps_atac_wt_rep_l[["W1_C3"]]
ht_list_WT_C <- hmaps_atac_wt_rep_l[["WT_C1"]] + hmaps_atac_wt_rep_l[["WT_C2"]] + hmaps_atac_wt_rep_l[["WT_C3"]]

ht_list_K1_V <- hmaps_atac_wt_rep_l[["K1_V1"]] + hmaps_atac_wt_rep_l[["K1_V2"]] + hmaps_atac_wt_rep_l[["K1_V3"]]
ht_list_KD_V <- hmaps_atac_wt_rep_l[["KD_V1"]] + hmaps_atac_wt_rep_l[["KD_V2"]] + hmaps_atac_wt_rep_l[["KD_V3"]]
ht_list_K1_C <- hmaps_atac_wt_rep_l[["K1_C1"]] + hmaps_atac_wt_rep_l[["K1_C2"]] + hmaps_atac_wt_rep_l[["K1_C3"]]
ht_list_KD_C <- hmaps_atac_wt_rep_l[["KD_C1"]] + hmaps_atac_wt_rep_l[["KD_C2"]] + hmaps_atac_wt_rep_l[["KD_C3"]]

#WT
png("../results/Human_cells/timepoint_reg/heatmaps_coverage/ht_list_wt_veh1h.png", width = 2, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + ht_list_W1_V, column_title = "           Control (1h)", 
     column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"), heatmap_legend_side = "bottom",
     heatmap_legend_list = wt_legend, gap = unit(1, "mm"), split = partition_wt)
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/ht_list_wt_veh36h.png", width = 2, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + ht_list_WT_V, column_title = "           Control (36h)", 
     column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"), heatmap_legend_side = "bottom",
     heatmap_legend_list = wt_legend, gap = unit(1, "mm"), split = partition_wt)
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/ht_list_wt_ccm1h.png", width = 2, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + ht_list_W1_C, column_title = "           Control CCM-like env (1h)", 
     column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"), heatmap_legend_side = "bottom",
     heatmap_legend_list = wt_legend, gap = unit(1, "mm"), split = partition_wt)
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/ht_list_wt_ccm36h.png", width = 2, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + ht_list_WT_C, column_title = "           Control CCM-like env (36h)", 
     column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"), heatmap_legend_side = "bottom",
     heatmap_legend_list = wt_legend, gap = unit(1, "mm"), split = partition_wt)
dev.off()

#KD
png("../results/Human_cells/timepoint_reg/heatmaps_coverage/ht_list_KD_veh1h.png", width = 2, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + ht_list_K1_V, column_title = "           siPDCD10 Veh (1h)", 
     column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"), heatmap_legend_side = "bottom",
     heatmap_legend_list = wt_legend, gap = unit(1, "mm"), split = partition_wt)
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/ht_list_KD_veh36h.png", width = 2, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + ht_list_KD_V, column_title = "           siPDCD10 Veh (36h)", 
     column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"), heatmap_legend_side = "bottom",
     heatmap_legend_list = wt_legend, gap = unit(1, "mm"), split = partition_wt)
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/ht_list_KD_ccm1h.png", width = 2, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + ht_list_K1_C, column_title = "           siPDCD10 CCM-like env (1h)", 
     column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"), heatmap_legend_side = "bottom",
     heatmap_legend_list = wt_legend, gap = unit(1, "mm"), split = partition_wt)
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/ht_list_KD_ccm36h.png", width = 2, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + ht_list_KD_C, column_title = "           siPDCD10 CCM-like env (36h)", 
     column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"), heatmap_legend_side = "bottom",
     heatmap_legend_list = wt_legend, gap = unit(1, "mm"), split = partition_wt)
dev.off()


## Junb bw ----

bw_file_chip <- list.files("../results/Human_cells/ChipSeq/bw/", full.names = T)[c(12, 11, 10, 9,
                                                                                   4, 3, 2, 1,
                                                                                   8, 7, 6, 5)]
bw_file_chip


bw_names_chip <- str_replace_all(basename(bw_file_chip), c("_sorted_rmdup.bw" = ""))

# names(bw_file) <- bw_names

bw_chip_list <- setNames(lapply(bw_file_chip, import), 
                         bw_names_chip)

# Normalize and save with correct filenames----
norm_chip_mtx_l <- normalize_bw_list(bw_chip_list, target = gr_sets_wt, extend = 1000)
# saveRDS(norm_chip_mtx_l, "../results/Human_cells/bw_norm_mtx/norm_chip_mtx_l.rds")
# norm_chip_mtx_l <- readRDS("../results/Human_cells/bw_norm_mtx/norm_chip_mtx_l.rds")

quantiles_junb_dep <- lapply(norm_chip_mtx_l[1:4], compute_quantiles)
quantiles_h3ac_dep <- lapply(norm_chip_mtx_l[5:8], compute_quantiles)
quantiles_h3me_dep <- lapply(norm_chip_mtx_l[9:12], compute_quantiles)

max(unlist(norm_chip_mtx_l[9:12]))

condition_names_chip <- c("Control (36h)", "Control\nCCM-like env (36h)", 
                         "siPDCD10\nVeh (36h)", "siPDCD10\nCCM-like env (36h)")

titles_hmaps_chip <- setNames(rep(condition_names_chip, 3),
                            names(norm_chip_mtx_l))


mtx_junb_wt_l <- map2(norm_chip_mtx_l[1:4], names(norm_chip_mtx_l)[1:4], 
                      ~{
                        .x %>%
                          as.data.frame() %>%
                          as.matrix()
                      }
)

mtx_h3ac_wt_l <- map2(norm_chip_mtx_l[5:8], names(norm_chip_mtx_l)[5:8], 
                      ~{
                        .x %>%
                          as.data.frame() %>%
                          as.matrix()
                      }
)

mtx_h3me_wt_l <- map2(norm_chip_mtx_l[9:12], names(norm_chip_mtx_l)[9:12], 
                      ~{
                        .x %>%
                          as.data.frame() %>%
                          as.matrix()
                      }
)


color_pal_map <- setNames(c(color_pal[1:2], color_pal[5:6]),
                          condition_names_chip)


annot_junb_l <- map2(titles_hmaps_chip, names(titles_hmaps_chip), 
                ~{
                  HeatmapAnnotation(Source = anno_block(gp = gpar(fill = color_pal_map[[.x]], 
                                                                  col = "black")), 
                                    height = unit(2.5, "mm"),
                                    show_legend = FALSE)
                }
)

col_fun1 <- circlize::colorRamp2(c(0, 70), c("white", "#66C2A5"))
col_fun3 <- circlize::colorRamp2(c(0, 60), c("white", "#CC075E"))
col_fun4 <- circlize::colorRamp2(c(0, 20), c("white", "dodgerblue"))

jun_legend <- Legend(col_fun = col_fun1, 
                    title = "RPKM", direction = "horizontal",
                    at = c(0,  35, 70),  # Define ticks on the legend
                    labels = c("0",  "35", "70"), title_position = "topcenter",
                    title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                    labels_gp = gpar(fontsize = 8, col = "black"),
                    legend_height = unit(2, "cm"),
                    legend_width = unit(2, "cm"),
                    title_gap = unit(1, "mm"))
hac_legend <- Legend(col_fun = col_fun3, 
                     title = "RPKM", direction = "horizontal",
                     at = c(0,  30, 60),  # Define ticks on the legend
                     labels = c("0",  "30", "60"), title_position = "topcenter",
                     title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                     labels_gp = gpar(fontsize = 8, col = "black"),
                     legend_height = unit(2, "cm"),
                     legend_width = unit(2, "cm"),
                     title_gap = unit(1, "mm"))
hme_legend <- Legend(col_fun = col_fun4, 
                     title = "RPKM", direction = "horizontal", title_position = "topcenter",
                     at = c(0,  10, 20),  # Define ticks on the legend
                     labels = c("0",  "10", "20"),
                     title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                     labels_gp = gpar(fontsize = 8, col = "black"),
                     legend_height = unit(2, "cm"),
                     legend_width = unit(2, "cm"),
                     title_gap = unit(1, "mm"))




## 🔹 Draw heatmaps ----
hmaps_jun_wt_l <- map2(mtx_junb_wt_l, names(mtx_junb_wt_l),
                        ~ {
                          annot <- annot_junb_l[[.y]]
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

hmaps_hac_wt_l <- map2(mtx_h3ac_wt_l, names(mtx_h3ac_wt_l),
                       ~ {
                         annot <- annot_junb_l[[.y]]
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

hmaps_hme_wt_l <- map2(mtx_h3me_wt_l, names(mtx_h3me_wt_l),
                       ~ {
                         annot <- annot_junb_l[[.y]]
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


ht_list_junbwt_wtregions_36h <-   hmaps_jun_wt_l[["junb_WT_VEH"]] +
  hmaps_jun_wt_l[["junb_WT_CCM"]] 

ht_list_junbkd_wtregions_36h <-   hmaps_jun_wt_l[["junb_KD_VEH"]] +
  hmaps_jun_wt_l[["junb_KD_CCM"]] 


ht_list_hk27wt_wtregions_36h <-   hmaps_hac_wt_l[["H3K27ac_WT_VEH"]] +
  hmaps_hac_wt_l[["H3K27ac_WT_CCM"]] 

ht_list_hk27kd_wtregions_36h <-   hmaps_hac_wt_l[["H3K27ac_KD_VEH"]] +
  hmaps_hac_wt_l[["H3K27ac_KD_CCM"]] 


ht_list_hkmewt_wtregions_36h <-   hmaps_hme_wt_l[["H3K4me1_WT_VEH"]] +
  hmaps_hme_wt_l[["H3K4me1_WT_CCM"]] 

ht_list_hkmekd_wtregions_36h <-   hmaps_hme_wt_l[["H3K4me1_KD_VEH"]] +
  hmaps_hme_wt_l[["H3K4me1_KD_CCM"]] 


## Save ATAC time points heatmap WT regions----

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_timepoints_wtjunb_wtregions36h.png", width = 1.25, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt2 + ht_list_junbwt_wtregions_36h, 
     #heatmap_legend_side = "bottom", heatmap_legend_list = combined_junb,
     gap = unit(1, "mm"), split = partition_wt,
     column_title = "      36 hours", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_timepoints_kdjunb_wtregions36h.png", width = 1.25, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt2 + ht_list_junbkd_wtregions_36h, 
     #heatmap_legend_side = "bottom", heatmap_legend_list = combined_junb,
     gap = unit(1, "mm"), split = partition_wt,
     column_title = "     36 hours", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()


png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_timepoints_wth3ac_wtregions36h.png", width = 1.25, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt2 + ht_list_hk27wt_wtregions_36h, 
     #heatmap_legend_side = "bottom", heatmap_legend_list = combined_junb,
     gap = unit(1, "mm"), split = partition_wt,
     column_title = "      36 hours", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_timepoints_kdh3ac_wtregions36h.png", width = 1.25, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt2 + ht_list_hk27kd_wtregions_36h, 
     #heatmap_legend_side = "bottom", heatmap_legend_list = combined_junb,
     gap = unit(1, "mm"), split = partition_wt,
     column_title = "     36 hours", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()


png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_timepoints_wth3me_wtregions36h.png", width = 1.25, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt2 + ht_list_hkmewt_wtregions_36h, 
     #heatmap_legend_side = "bottom", heatmap_legend_list = combined_junb,
     gap = unit(1, "mm"), split = partition_wt,
     column_title = "     36 hours", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_timepoints_kdh3me_wtregions36h.png", width = 1.25, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt2 + ht_list_hkmekd_wtregions_36h, 
     #heatmap_legend_side = "bottom", heatmap_legend_list = combined_junb,
     gap = unit(1, "mm"), split = partition_wt,
     column_title = "     36 hours", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()


### hmaps atac counts----

gr_sets_wt

cts <- read_csv("../results/Human_cells/D3_ATAC/data/log2cpm.csv", col_select = -1)

cts_filt <- cts %>% 
  select(peaks,
         starts_with("W1_V"), starts_with("W1_C"), 
         starts_with("WT_V"), starts_with("WT_C"), 
         starts_with("K1_V"), starts_with("K1_c"), 
         starts_with("KD_V"), starts_with("KD_C")) %>% 
  filter(peaks %in% gr_sets_wt$peaks) %>% 
  column_to_rownames("peaks")

mtx_scaled <-  t(scale(t(cts_filt)))

col_fun_atac <- circlize::colorRamp2(seq(-3, 3, length = 100), viridis(100))
cts_legend <- Legend(col_fun = col_fun_atac, direction = "horizontal",
                     title = "Z-score log2 (CPM +1)",  title_position = "topcenter",
                     title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                     labels_gp = gpar(fontsize = 8, col = "black"),
                     legend_height = unit(2, "cm"),
                     title_gap = unit(1, "mm"))

combined_legend_cts <- packLegend(cts_legend, direction = "vertical")

heatmap_builder <- function(mtx, order_rows, name = NULL, annotation){
  
  ComplexHeatmap::Heatmap(mtx, 
                          col = col_fun_atac,
                          show_row_names = FALSE,
                          show_heatmap_legend = FALSE, 
                          column_title = name, 
                          column_title_gp = gpar(fontsize = 7),
                          show_column_names = FALSE,
                          cluster_rows = FALSE, border = TRUE,
                          top_annotation = annotation,
                          cluster_columns = FALSE, 
                          row_order = order_rows)
}

wt_1hv <- heatmap_builder(mtx_scaled[ gr_sets_wt$peaks, c(1:3)], order_for_heatmaps, annotation = annot_l[["W1_V1"]])
wt_1hc <- heatmap_builder(mtx_scaled[ gr_sets_wt$peaks, c(4:6)], order_for_heatmaps, annotation = annot_l[["W1_C1"]])
wt_36hv <- heatmap_builder(mtx_scaled[ gr_sets_wt$peaks, c(7:9)], order_for_heatmaps, annotation = annot_l[["WT_V1"]])
wt_36hc <- heatmap_builder(mtx_scaled[ gr_sets_wt$peaks, c(10:12)], order_for_heatmaps, annotation = annot_l[["WT_C1"]])

kd_1hv <- heatmap_builder(mtx_scaled[ gr_sets_wt$peaks, c(13:15)], order_for_heatmaps, annotation = annot_l[["K1_V1"]])
kd_1hc <- heatmap_builder(mtx_scaled[ gr_sets_wt$peaks, c(16:18)], order_for_heatmaps, annotation = annot_l[["K1_C1"]])
kd_36hv <- heatmap_builder(mtx_scaled[ gr_sets_wt$peaks, c(19:21)], order_for_heatmaps, annotation = annot_l[["KD_V1"]])
kd_36hc <- heatmap_builder(mtx_scaled[ gr_sets_wt$peaks, c(22:24)], order_for_heatmaps, annotation = annot_l[["KD_C1"]])


cts_wt_wtregions_1h <- wt_1hv +
  wt_1hc
cts_wt_wtregions_36h <- wt_36hv +
  wt_36hc 

cts_kd_wtregions_1h <- kd_1hv +
  kd_1hc
cts_kd_wtregions_36h <- kd_36hv +
  kd_36hc 

## Save ATAC time points heatmap cts----
png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_cts_wt_wtregions_1h.png", width = 1.45, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + cts_wt_wtregions_1h,  gap = unit(1, "mm"), show_heatmap_legend = FALSE, 
     heatmap_legend_list = combined_legend_cts, split = partition_wt, row_title_gp = gpar(fontsize = 9),
     column_title = "         1 hour", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_cts_wt_wtregions_36h.png", width = 1.25, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt2 + cts_wt_wtregions_36h,  gap = unit(1, "mm"), show_heatmap_legend = FALSE, 
     heatmap_legend_list = combined_legend_cts, split = partition_wt, row_title_gp = gpar(fontsize = 9),
     column_title = "     36 hours", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()


png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_cts_kd_wtregions_1h.png", width = 1.45, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt + cts_kd_wtregions_1h,  gap = unit(1, "mm"), show_heatmap_legend = FALSE, 
     heatmap_legend_list = combined_legend_cts, split = partition_wt, row_title_gp = gpar(fontsize = 9),
     column_title = "         1 hour", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()
png("../results/Human_cells/timepoint_reg/heatmaps_coverage/atac_cts_kd_wtregions_36h.png", width = 1.25, height = 3.14, units = "in", res = 300)
draw(partition_hp_wt2 + cts_kd_wtregions_36h,  gap = unit(1, "mm"), show_heatmap_legend = FALSE, 
     heatmap_legend_list = combined_legend_cts, split = partition_wt, row_title_gp = gpar(fontsize = 9),
     column_title = "     36 hours", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()


combined_legend_wt_cts <- packLegend(cts_legend,  direction = "horizontal")
combined_legend_kd_cts <- packLegend(cts_legend,  direction = "horizontal")

png(filename = "../results/Human_cells/timepoint_reg/heatmaps_coverage/legendwt_cts.png", width = 4, height = 1, units = "in", res = 300)
ComplexHeatmap::draw(combined_legend_wt_cts)
dev.off()
png(filename = "../results/Human_cells/timepoint_reg/heatmaps_coverage/legendkd_cts.png", width = 4, height = 1, units = "in", res = 300)
ComplexHeatmap::draw(combined_legend_kd_cts)
dev.off()

### hmaps RNA counts----

cts_rna <- read_csv("../results/Human_cells/D3_RNA/data/log2cpm_filt.csv")
gene_ehn_pairs <- read_csv("../results/Human_cells/D3_ATAC/data/pcc_ccre_genes_500kb_filt.csv")

regions_for_RNA <- data.frame(peaks = c(early_only_WT_strict$peaks, persistent_WT_strict$peaks, late_only_WT_strict$peaks),
                              Classification = fct_inorder(rep(c("Early response",
                                                                 #"Early response",
                                                                 "Late response"), 
                                                               c(length(early_only_WT_strict$peaks) + length(persistent_WT_strict$peaks),
                                                                 length(late_only_WT_strict$peaks)))))

cts_rna_filt <- cts_rna %>% 
  select(gene,
         starts_with("WT_V"), starts_with("WT_C"), 
         starts_with("KD_V"), starts_with("KD_C")) %>% 
  column_to_rownames("gene")

cts_rna_filt_scaled <-  as.data.frame(t(scale(t(cts_rna_filt)))) %>% 
  rownames_to_column("gene") %>%
  drop_na()


cts_enh_RNA <- regions_for_RNA %>% 
  inner_join(gene_ehn_pairs %>% select(peaks, gene), by = "peaks") %>% 
  inner_join(cts_rna_filt_scaled, by = "gene")

write_csv(cts_enh_RNA, "../results/Human_cells/timepoint_reg/tables/cts_enh_RNA_D3.csv")

max(cts_enh_RNA[-c(1:3)])
min(cts_enh_RNA[-c(1:3)])

col_fun_rna <- circlize::colorRamp2(seq(-3, 3, length = 100), plasma(100))
cts_rna_legend <- Legend(col_fun = col_fun_rna, direction = "horizontal",
                     title = "Z-score log2 (CPM +1)", title_position = "topcenter",
                     title_gp = gpar(fontsize = 8.5, fontface = "bold"),   # Legend title appearance
                     labels_gp = gpar(fontsize = 8, col = "black"),
                     legend_height = unit(2, "cm"),
                     title_gap = unit(1, "mm"))

combined_legend_rna_cts <- packLegend(cts_rna_legend, direction = "vertical")
partition_rna <- fct_inorder(rep(c("Early response",
                                   #"Early response",
                                   "Late response"), 
                                 c(length(cts_enh_RNA$Classification[cts_enh_RNA$Classification == "Early response"]), 
                                  # length(cts_enh_RNA$Classification[cts_enh_RNA$Classification == "Early response"]),
                                   length(cts_enh_RNA$Classification[cts_enh_RNA$Classification == "Late response"]))))

partition_rna_wt <- Heatmap(partition_rna, col=structure(c("#EABFCB", 
                                                           #"#C191A1", 
                                                           "#A4508B"), names = c("Early response",
                                                                                 #"Early response",
                                                                                 "Late response")), 
                           show_column_names = FALSE, name = " ", row_title = c("",""),
                           show_heatmap_legend = FALSE, row_title_gp = gpar(fontsize = 0), 
                           show_row_names = FALSE, width=unit(1,'mm'))



# Right-side annotation
ha <- rowAnnotation(motif_labels = anno_mark(at = c(36, 69, 124, 307, 679, 793),
                                             labels = c("CCL2", "ICAM1", "SERPINE1", "CCL2", "VCAM1", "VEGFA"),
                                             labels_gp = gpar(fontsize = 8, col = "black"),
                                             padding = unit(1, "mm")))

heatmap_builder_rna <- function(mtx, order_rows, name = NULL, annotation, right_annot = NULL){
  
  ComplexHeatmap::Heatmap(mtx, 
                          col = col_fun_rna,
                          show_row_names = FALSE,
                          show_heatmap_legend = FALSE, 
                          column_title = name, 
                          column_title_gp = gpar(fontsize = 7),
                          show_column_names = FALSE,
                          cluster_rows = FALSE, right_annotation = right_annot,
                          top_annotation = annotation,
                          cluster_columns = FALSE, 
                          row_order = order_rows)
}


wt_36hv_rna <- heatmap_builder_rna(cts_enh_RNA[ , c(4:6)], NULL, annotation = annot_l[["WT_V1"]])
wt_36hc_rna <- heatmap_builder_rna(cts_enh_RNA[ , c(7:9)], NULL, annotation = annot_l[["WT_C1"]], right_annot = ha)

kd_36hv_rna <- heatmap_builder_rna(cts_enh_RNA[ , c(10:12)], NULL, annotation = annot_l[["KD_V1"]])
kd_36hc_rna <- heatmap_builder_rna(cts_enh_RNA[ , c(13:15)], NULL, annotation = annot_l[["KD_C1"]], right_annot = ha)


cts_wt_wtregions_36h <- wt_36hv_rna +
  wt_36hc_rna 

cts_kd_wtregions_36h <- kd_36hv_rna +
  kd_36hc_rna 
## Save RNA time points heatmap cts----

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/rna_cts_wt_wtregions_36h.png", width = 2.05, height = 3.14, units = "in", res = 300)
draw(partition_rna_wt + cts_wt_wtregions_36h,  gap = unit(1, "mm"), show_heatmap_legend = FALSE, 
     heatmap_legend_list = combined_legend_rna_cts, split = partition_rna, row_title_gp = gpar(fontsize = 9),
     column_title = "36 hours             ", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/rna_cts_kd_wtregions_36h.png", width = 2.05, height = 3.14, units = "in", res = 300)
draw(partition_rna_wt + cts_kd_wtregions_36h,  gap = unit(1, "mm"), show_heatmap_legend = FALSE, 
     heatmap_legend_list = combined_legend_rna_cts, split = partition_rna, row_title_gp = gpar(fontsize = 9),
     column_title = "36 hours             ", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()


combined_legend <- packLegend(cts_legend, jun_legend, hac_legend, hme_legend, cts_rna_legend, 
                              direction = "horizontal", gap = unit(0.45, "in"))

png(filename = "../results/Human_cells/timepoint_reg/heatmaps_coverage/legend_wt.png", width = 6.8, height = 0.5, units = "in", res = 300)
ComplexHeatmap::draw(combined_legend, x = unit(0.5, "npc"))
dev.off()


# Optionally save session info
sink("../logs/sessioninfo_human_16_coverage_D3_atac_h3k27ac_me1.txt")
sessionInfo()
sink()



