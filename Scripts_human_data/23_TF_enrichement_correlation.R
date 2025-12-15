pacman::p_load("dplyr", "tidyr", "purrr", "ggrepel", "ggpubr",
               "ggplot2", "stringr", "ggtext",
               "patchwork", "ggpubr", "glue",
               "vroom", "ggh4x", "colorspace")

plot_folder <- "../results/human_mouse_integration/plots/TF_corr"
dir.create(plot_folder, showWarnings = FALSE, recursive = TRUE)
motif_enrichment_per_cell <- read.csv("../results/RNA_ATAC_integration/Homer/results_tables/enhancer_ko_cCRE_knownTF_percell.csv")

motif_enrichment_CCM3KD <- read.csv("../results/Human_cells/Homer/results_tables/TF_enrichment_up_FDR_05.csv")

cell_types <- unique(motif_enrichment_per_cell$cell_type)
names_plot <- c("Artery <i>PDCD10</i><sup>BECKO</sup>", 
               "Capillary artery <i>PDCD10</i><sup>BECKO</sup>",
               "Capillary <i>PDCD10</i><sup>BECKO</sup>",
               "Capillary vein <i>PDCD10</i><sup>BECKO</sup>", 
               "Large vein <i>PDCD10</i><sup>BECKO</sup>")
names(names_plot) <- cell_types

# Create a named list of plots
plot_list <- lapply(cell_types, function(ct) {
  df_ko <- motif_enrichment_per_cell %>%
    filter(cell_type == ct) %>%
    rename(log10_pvalue_KO = log10_pvalue)
  
  df_kd <- motif_enrichment_CCM3KD %>%
    filter(condition == "KD_CCM_vs_WT_VEH") %>%
    dplyr::rename(log10_pvalue_KD_CCM = log10_pvalue)
  
  merged_df <- inner_join(df_ko, df_kd, by = "TF_motif")
  
  top_motifs <-  merged_df %>%
    filter(TF_motif %in% c("JunB", "ERG", "Fosl2", "Etv2", "Fos", "ETS1"))
  
  name_plot <- names_plot[[ct]]
  
  p <- ggplot(merged_df, aes(x = log10_pvalue_KO, y = log10_pvalue_KD_CCM)) +
    geom_point(color = "black", size = 0.7) +
    geom_point(data = top_motifs, aes(x = log10_pvalue_KO, y = log10_pvalue_KD_CCM),
               color = "darkmagenta", size = 0.7) +
    geom_text_repel(data = top_motifs, aes(label = TF_motif), fontface = "bold", size = 3.5, 
                    color = "darkmagenta", min.segment.length = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue", size = 0.3, linetype = "dashed") +
    stat_cor(method = "pearson",
             label.x = max(merged_df$log10_pvalue_KO) * 0.6,
             label.y = max(merged_df$log10_pvalue_KD_CCM) * 0.1,
             size = 3,
             fontface = "bold") +
    labs(title = name_plot, 
         x = glue("-log10 p-value<br>{name_plot} (mouse)"),
         y = "-log10 p-value<br>siPDCD10 + CCM-like env (human)") +
    scale_x_continuous(limits = c(0, 330)) +
    scale_y_continuous(limits = c(0, 330)) +
    theme_bw(base_size = 10) + theme(title = element_markdown(), axis.title.x = element_markdown())
  
  p_save <- ggplot(merged_df, aes(x = log10_pvalue_KO, y = log10_pvalue_KD_CCM)) +
    geom_point(color = "black", size = 0.7) +
    geom_point(data = top_motifs, aes(x = log10_pvalue_KO, y = log10_pvalue_KD_CCM),
               color = "darkmagenta", size = 0.7) +
    geom_text_repel(data = top_motifs, aes(label = TF_motif), fontface = "bold", size = 3.5, 
                    color = "darkmagenta", min.segment.length = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue", size = 0.3, linetype = "dashed") +
    stat_cor(method = "pearson",
             label.x = max(merged_df$log10_pvalue_KO) * 0.6,
             label.y = max(merged_df$log10_pvalue_KD_CCM) * 0.1,
             size = 3,
             fontface = "bold") +
    labs(title = "TF motif enrichment", 
         x = glue("-log10 p-value<br>{name_plot} (mouse)"),
         y = "-log10 p-value<br>siPDCD10 + CCM-like env (human)") +
    scale_x_continuous(limits = c(0, max(merged_df$log10_pvalue_KO) * 1.1)) +
    scale_y_continuous(limits = c(0, max(merged_df$log10_pvalue_KD_CCM) * 1.1)) +
    theme_bw() + 
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_markdown(face = "bold"))

  ggsave(glue("{plot_folder}/TF_corr_{ct}.png"),
         plot = p_save, width = 3.2, height = 3.2, unit = "in", device = "png")

  return(p)
})

# Combine all plots
plots_wrapped <- wrap_plots(plot_list, ncol = 2) + plot_annotation(title = "TF motif correlation across vascular subtypes")


ggsave(glue("{plot_folder}/TF_corr_across_cell_types.png"), 
       plot = plots_wrapped, width = 10, height = 7, unit = "in", device = "png")

# Optionally save session info
sink("../logs/sessioninfo_human_23_TF_enrichement_correlation.txt")
sessionInfo()
sink()


