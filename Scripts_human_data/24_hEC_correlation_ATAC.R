
pacman::p_load("dplyr", "tidyr", "ggpubr", "glue", "stringr",
               "ggtext", "tibble", "hexbin", "MASS", "readr")

plot_dir <- "../results/human_mouse_integration/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE)

cts_d3_p <- read.csv("../results/Human_cells/D3_ATAC/data/log2cpm.csv", row.names = 1) %>% 
  dplyr::select(peaks, contains(("WT")))

cts_ae_p <- read.csv("../results/Human_cells/HPAEC_ATAC/data/log2cpm.csv")

cts_ve_p <- read.csv("../results/Human_cells/HUVEC_ATAC/data/log2cpm.csv")


cts_d3_gr <- cts_d3_p %>% 
  column_to_rownames("peaks") %>%
  dplyr::select(matches("_C1|_C2|_C3")) %>%
  mutate(avg_cpm_hu = rowMeans(.)) %>%
  rownames_to_column("to_split") %>% 
  mutate(peaks = to_split) %>% 
  separate(to_split, into = c("chr", "start", "end"), 
           sep = "-", remove = TRUE ) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

cts_ae_gr <- cts_ae_p %>% 
  column_to_rownames("peaks") %>%
  dplyr::select(matches("_C1|_C2|_C3")) %>%
  mutate(avg_cpm_hu = rowMeans(.)) %>%
  rownames_to_column("to_split") %>% 
  mutate(peaks = to_split) %>%
  separate(to_split, into = c("chr", "start", "end"), 
           sep = "-", remove = TRUE) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

cts_ve_gr <- cts_ve_p %>% 
  column_to_rownames("peaks") %>%
  dplyr::select(matches("_C1|_C2|_C3")) %>%
  mutate(avg_cpm_hu = rowMeans(.)) %>%
  rownames_to_column("to_split") %>% 
  mutate(peaks = to_split) %>%
  separate(to_split, into = c("chr", "start", "end"), 
           sep = "-", remove = TRUE) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

hits_d3_hpaec <- findOverlaps(cts_d3_gr, cts_ae_gr)
hits_d3_huvec <- findOverlaps(cts_d3_gr, cts_ve_gr)


common_d3_ae_df <- data.frame(hCMEC_set1 = cts_d3_gr[queryHits(hits_d3_hpaec)]$peaks,
                              HPAEC_set2 = cts_ae_gr[subjectHits(hits_d3_hpaec)]$peaks,
                              avgCPM_hCMEC = cts_d3_gr[queryHits(hits_d3_hpaec)]$avg_cpm_hu,
                              avgCPM_HPAEC = cts_ae_gr[subjectHits(hits_d3_hpaec)]$avg_cpm_hu)

common_d3_ve_df <- data.frame(hCMEC_set1 = cts_d3_gr[queryHits(hits_d3_huvec)]$peaks,
                              HUVEC_set2 = cts_ve_gr[subjectHits(hits_d3_huvec)]$peaks,
                              avgCPM_hCMEC = cts_d3_gr[queryHits(hits_d3_huvec)]$avg_cpm_hu,
                              avgCPM_HUVEC = cts_ve_gr[subjectHits(hits_d3_huvec)]$avg_cpm_hu)


dens_d3_ae <- kde2d(common_d3_ae_df$avgCPM_hCMEC, common_d3_ae_df$avgCPM_HPAEC, n = 200)
common_d3_ae_df$density <- fields::interp.surface(dens_d3_ae, cbind(common_d3_ae_df$avgCPM_hCMEC, common_d3_ae_df$avgCPM_HPAEC))

dens_d3_ve <- kde2d(common_d3_ve_df$avgCPM_hCMEC, common_d3_ve_df$avgCPM_HUVEC, n = 200)
common_d3_ve_df$density <- fields::interp.surface(dens_d3_ve, cbind(common_d3_ve_df$avgCPM_hCMEC, common_d3_ve_df$avgCPM_HUVEC))

# Plot
hCMEC_HPAEC <- ggplot(common_d3_ae_df, aes(x = avgCPM_hCMEC, y = avgCPM_HPAEC)) +
  geom_point(aes(color = density), size = 0.2, alpha = 0.7) +
  scale_color_viridis_c(option = "C") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  labs(title = "Chromatin accessibility correlation",
       x = "log2(CPM+1)<br> hCMEC + TNF/DMOG",
       y = "log2(CPM+1)<br> HPAEC + TNF/DMOG",
       color = "Density") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_markdown(face = "bold", lineheight = 1.2)) +
  stat_cor(method = "pearson", label.x = 2, label.y = max(common_d3_ae_df$avgCPM_HPAEC), size = 4)


# Plot
hCMEC_HUVEC <- ggplot(common_d3_ve_df, aes(x = avgCPM_hCMEC, y = avgCPM_HUVEC)) +
  geom_point(aes(color = density), size = 0.2, alpha = 0.7) +
  scale_color_viridis_c(option = "C") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  labs(title = "Chromatin accessibility correlation",
       x = "log2(CPM+1)<br> hCMEC + TNF/DMOG",
       y = "log2(CPM+1)<br> HUVEC + TNF/DMOG",
       color = "Density") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_markdown(face = "bold", lineheight = 1.2)) +
  stat_cor(method = "pearson", label.x = 2, label.y = max(common_d3_ve_df$avgCPM_HUVEC), size = 4)

ggsave(glue("{plot_dir}/hCMEC_HPAEC_peak_corr.png"),
       plot = hCMEC_HPAEC, width = 3.5, height = 3, dpi = 300)

ggsave(glue("{plot_dir}/hCMEC_HUVEC_peak_corr.png"),
       plot = hCMEC_HUVEC, width = 3.5, height = 3, dpi = 300)

# Optionally save session info
sink("../logs/sessioninfo_human_24_hEC_correlation_ATAC.txt")
sessionInfo()
sink()
