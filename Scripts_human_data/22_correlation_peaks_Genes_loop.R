
pacman::p_load("dplyr", "tidyr", "ggpubr", "glue", "stringr",
               "ggtext", "tibble", "hexbin", "MASS")

plot_dir <- "../results/human_mouse_integration/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE)

cts_hu_g <- read.csv("../results/Human_cells/D3_RNA/data/log2cpm_filt.csv")
cts_hu_p <- read.csv("../results/Human_cells/D3_ATAC/data/log2cpm.csv", row.names = 1)

cts_ms_g <- read.csv("../results/RNA/data/log2cpm_filtered.csv")
cts_ms_p <- read.csv("../results/ATAC/data/log2cpm.csv")


mouse_human_peaks <- read.csv("../results/human_mouse_integration/data/hu_ms_ortholog_peaks.csv") 

jax_ortholog_path <- "utils/jax_ms_hu_orthologs.tsv"
if (!file.exists(jax_ortholog_path)) {
  download.file("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
                destfile = jax_ortholog_path)
}
mouse_human_genes <- read.table(jax_ortholog_path, sep = "\t", header = TRUE)

convert_mouse_to_human <- function(gene_list){
  mouse_filtered <- mouse_human_genes %>% filter(Symbol %in% gene_list, Common.Organism.Name == "mouse, laboratory")
  human_filtered <- mouse_human_genes %>% filter(Common.Organism.Name == "human")
  merged <- left_join(mouse_filtered, human_filtered, by = "DB.Class.Key") %>%
    dplyr::select(Symbol.x, Symbol.y) %>%
    dplyr::rename(gene_ms = Symbol.x, gene_hu = Symbol.y)
  return(merged)
}

subtypes <- c("Artery", "CapArt", "Cap", "CapVein", "LargeVein")
KO_names <- gsub("\\d_.*$", "", colnames(cts_ms_g)[grep("_KO", colnames(cts_ms_g))]) %>% unique()
names(subtypes) <- KO_names

# Human column suffixes
human_rna_cols <- c("KD_C1", "KD_C2", "KD_C3")
human_atac_cols <- c("KD_C1", "KD_C2", "KD_C3")

conserved_per_cell_type <- list()
for (sub in names(subtypes)) {
  message(glue("▶ Processing {sub}"))
  
  file_suffix <- subtypes[[sub]]
  cell_type <- str_to_sentence(str_replace_all(sub, c("_KO" = "", "_" = " ")))
    
    gsub("_KO", "", sub)
  ## --- Load DEG and DAR
  deg_path <- glue("../results/RNA/edgeR/deg_difexp_tables/{file_suffix}_KO_vs_Flox_diffexp.csv")
  dar_path <- glue("../results/ATAC/edgeR/dar_diffacc_tables/{file_suffix}_KO_vs_Flox_diffacc.csv")
  deg <- read.csv(deg_path)
  dar <- read.csv(dar_path)
  
  ## --- Gene ortholog conversion
  mouse_to_human <- convert_mouse_to_human(deg$gene)
  
  ## --- RNA average CPM
  rna_cols <- grep(glue("^{sub}"), colnames(cts_ms_g), value = TRUE)
  ms_rna_avg <- cts_ms_g %>%
    column_to_rownames("gene") %>%
    dplyr::select(all_of(rna_cols)) %>%
    mutate(avg_cpm_ms = rowMeans(.)) %>%
    rownames_to_column("gene_ms")
  
  hu_rna_avg <- cts_hu_g %>%
    column_to_rownames("gene") %>%
    dplyr::select(all_of(human_rna_cols)) %>%
    mutate(avg_cpm_hu = rowMeans(.)) %>%
    rownames_to_column("gene_hu")
  
  gene_df <- deg %>%
    rename(gene_ms = gene) %>%
    inner_join(mouse_to_human, by = "gene_ms") %>%
    inner_join(ms_rna_avg, by = "gene_ms") %>%
    inner_join(hu_rna_avg, by = "gene_hu")
  
  # Compute 2D density
  dens <- kde2d(gene_df$avg_cpm_hu, gene_df$avg_cpm_ms, n = 200)
  gene_df$density <- fields::interp.surface(dens, cbind(gene_df$avg_cpm_hu, gene_df$avg_cpm_ms))
  
  # Plot
  p_gene <- ggplot(gene_df, aes(x = avg_cpm_hu, y = avg_cpm_ms)) +
    geom_point(aes(color = density), size = 0.2, alpha = 0.7) +
    scale_color_viridis_c(option = "C") +
    geom_smooth(method = "lm", color = "red", linetype = "dashed") +
    labs(title = glue("Gene correlation"),
         x = glue("log2(CPM+1)<br>{cell_type} <i>PDCD10</i><sup>BECKO</sup> (mouse)"),
         y = "log2(CPM+1)<br>siPDCD10 + CCM-like env (human)",
         color = "Density") +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_markdown(face = "bold")) +
    stat_cor(method = "pearson", label.x = 2, label.y = max(gene_df$avg_cpm_ms), size = 4)
  
  # ggsave(glue("{plot_dir}/{sub}_gene_corr.png"),
  #        plot = p_gene, width = 3.5, height = 3.5, dpi = 300)
  
  ## --- Peaks
  peak_cols <- grep(glue::glue("^{sub}"), colnames(cts_ms_p), value = TRUE)
  ms_peak_avg <- cts_ms_p %>%
    column_to_rownames("peaks") %>%
    dplyr::select(all_of(peak_cols)) %>%
    mutate(avg_cpm_ms = rowMeans(.)) %>%
    rownames_to_column("peaks_ms")
  
  hu_peak_avg <- cts_hu_p %>%
    column_to_rownames("peaks") %>%
    dplyr::select(all_of(human_atac_cols)) %>%
    mutate(avg_cpm_hu = rowMeans(.)) %>%
    rownames_to_column("peaks_hu")
  
  peak_df <- dar %>%
    rename(peaks_ms = peaks) %>%
    inner_join(mouse_human_peaks, by = "peaks_ms") %>%
    inner_join(ms_peak_avg, by = "peaks_ms") %>%
    inner_join(hu_peak_avg, by = "peaks_hu")
  
  # Compute 2D density
  dens <- kde2d(peak_df$avg_cpm_hu, peak_df$avg_cpm_ms, n = 200)
  peak_df$density <- fields::interp.surface(dens, cbind(peak_df$avg_cpm_hu, peak_df$avg_cpm_ms))
  
  # Plot
  p_peak <- ggplot(peak_df, aes(x = avg_cpm_hu, y = avg_cpm_ms)) +
    geom_point(aes(color = density), size = 0.2, alpha = 0.7) +
    scale_color_viridis_c(option = "C") +
    geom_smooth(method = "lm", color = "red", linetype = "dashed") +
    labs(title = glue("Chromatin accessibility correlation"),
         x = glue("log2(CPM+1)<br>{cell_type} <i>PDCD10</i><sup>BECKO</sup> (mouse)"),
         y = "log2(CPM+1)<br>siPDCD10 + CCM-like env (human)",
         color = "Density") +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_markdown(face = "bold")) +
    stat_cor(method = "pearson", label.x = 2, label.y = max(peak_df$avg_cpm_ms), size = 4)
  
  # ggsave(glue("{plot_dir}/{sub}_peak_corr.png"),
  #        plot = p_peak, width = 3.5, height = 3.5, dpi = 300)
  # 
  conserved_per_cell_type[[sub]] <- peak_df %>% dplyr::select(peaks_hu, starts_with(c("logFC", "FDR", "diff")))
}


peaks_conserved_changing <- reduce(conserved_per_cell_type, full_join, by = "peaks_hu") 

write.csv(peaks_conserved_changing, "../results/human_mouse_integration/data/hu_ms_conserved_changing.csv", row.names = FALSE)



# Optionally save session info
sink("../logs/sessioninfo_human_22_correlation_peaks_Genes_loop_ms_d3.txt")
sessionInfo()
sink()
















