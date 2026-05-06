# Script 09: Differential accessibility with edgeR from pseudobulk ATAC
pacman::p_load("edgeR", "limma", "ggplot2", "dplyr", "glue", "readr", "stringr",
               "tidyr", "tibble", "purrr", "EnhancedVolcano", "ComplexUpset", "ggtext")

# Output folders----
dar_folder <- "../results/ATAC/edgeR/dar_tables"
dar_difexp_folder <- "../results/ATAC/edgeR/dar_diffacc_tables"
rds_folder <- "../results/ATAC/rds"
volcano_dir <- "../results/plots/graphs_ATAC"
summary_dir <- "../results/plots/graphs_endo"

dir.create(dar_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(dar_difexp_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_folder, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# Directories
cts_file <- "../results/ATAC/data/counts_atac_all.csv"
peakset_flox_file <- "../results/ATAC/filtered_peaksets/peakset_all_flox/all_flox_filtered_union.final.peakset"
peakset_ko_dir <- "../results/ATAC/filtered_peaksets/peakset_by_celltype"

# Load counts
cts <- read.csv(cts_file, row.names = 1)
colnames(cts) <- gsub("\\.", "_", colnames(cts))
colnames(cts) <- paste0(colnames(cts), "_ATAC")

# Calculate CPM and log2CPM
cpm <- edgeR::cpm(cts)
log2cpm <- log2(cpm + 1)
# Save CPM tables
write_csv(rownames_to_column(as.data.frame(cpm), var = "peaks"), "../results/ATAC/data/cpm.csv")
write_csv(rownames_to_column(as.data.frame(log2cpm), var = "peaks"), "../results/ATAC/data/log2cpm.csv")

# --- Flox-only analysis (subset & master peak filtering)
flox_cols <- colnames(cts)[str_detect(colnames(cts), "Flox")]
cts_flox <- cts[ , flox_cols]
flox_peaks <- read.table(peakset_flox_file, header = TRUE) %>%
  mutate(PeakID = paste(chr, start, end, sep = "-")) %>%
  pull(PeakID)


groups_flox <- c(rep("AF", 3), rep("cAF", 3), rep("cF", 3), rep("cVF", 3), rep("LVF", 3))
dge_flox <- DGEList(counts = cts_flox, group = factor(groups_flox))
dge_flox <- dge_flox[flox_peaks[flox_peaks %in% rownames(dge_flox)], ]
dge_flox$samples$lib.size <- colSums(dge_flox$counts)
dge_flox <- calcNormFactors(dge_flox)

#Optional MDS
# plotMDS(dge_flox, method = "bcv", col = as.numeric(dge_flox$samples$group), pch = 20)
# legend("topleft", legend = levels(dge_flox$samples$group), col = 1:5, pch = 20)

design_flox <- model.matrix(~ 0 + dge_flox$samples$group)
colnames(design_flox) <- levels(dge_flox$samples$group)
dge_flox <- estimateDisp(dge_flox, design_flox)

# Optional BCV plot
# plotBCV(dge_flox)

fit_flox <- glmFit(dge_flox, design_flox)
contrasts_flox <- makeContrasts(Artery_vs_rest    = AF - (cAF + cVF + cF + LVF)/4,
                                CapArt_vs_rest    = cAF - (AF + cVF + cF + LVF)/4,
                                Cap_vs_rest       = cF - (AF + cAF + cVF + LVF)/4,
                                CapVein_vs_rest   = cVF - (AF + cAF + cF + LVF)/4,
                                LargeVein_vs_rest = LVF - (AF + cAF + cVF + cF)/4,
                                levels = design_flox)

# Run DE for Flox
dar.l <- list()
for (i in colnames(contrasts_flox)) {
  lrt <- glmLRT(fit_flox, contrast = contrasts_flox[, i])
  tab <- topTags(lrt, n = Inf)$table %>%
    rownames_to_column("peaks") %>%
    rename_with(~ paste0(., "_", i), -peaks)
  dar.l[[i]] <- tab
  write_csv(tab, file.path(dar_folder, paste0(i, ".csv")))
}

# KO differential accessibility analysis
group_defs <- list(artery             = rep(c("AF", "AK"), each = 3),
                   capillary_artery = rep(c("cAF", "cAK"), each = 3),
                   capillary          = rep(c("cF", "cK"), each = 3),
                   capillary_vein   = rep(c("cVF", "cVK"), each = 3),
                   large_vein         = rep(c("LVF", "LVK"), each = 3))

contrasts_ko <- c(Artery_KO_vs_Flox        = "AK - AF",
                  CapArt_KO_vs_Flox        = "cAK - cAF",
                  Cap_KO_vs_Flox           = "cK - cF",
                  CapVein_KO_vs_Flox       = "cVK - cVF",
                  LargeVein_KO_vs_Flox     = "LVK - LVF")

contrast_numbers <- list(artery = 1,
                         capillary_artery = 2,
                         capillary = 3,
                         capillary_vein = 4,
                         large_vein = 5)


for (celltype in names(group_defs)) {
  peaks_file <- file.path(peakset_ko_dir, paste0(celltype, ".filteredNfixed.union.final.peakset"))
  peaks <- read.table(peaks_file, header = TRUE) %>%
    mutate(PeakID = paste(chr, start, end, sep = "-")) %>%
    pull(PeakID)
  
  cts_subset <- cts %>%
    select(matches(paste0("^", celltype, "_(Flox|KO)[0-9]+_ATAC$")))
  cts_subset <- cts_subset[rownames(cts_subset) %in% peaks, ]
  
  d <- DGEList(counts = cts_subset, group = factor(group_defs[[celltype]]))
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  design <- model.matrix(~ 0 + d$samples$group)
  colnames(design) <- levels(d$samples$group)
  
  d <- estimateDisp(d, design)
  fit <- glmFit(d, design)
  
  i <- contrast_numbers[[celltype]]
  contrast_mat <- makeContrasts(contrasts = contrasts_ko[i], levels = design)
  contrast_name <- names(contrasts_ko[i])
  
  lrt <- glmLRT(fit, contrast = contrast_mat)
  tab <- topTags(lrt, n = Inf)$table %>%
    rownames_to_column("peaks") %>%
    rename_with(~ paste0(., "_", contrast_name), -peaks)
  dar.l[[contrast_name]] <- tab
  write_csv(tab, file.path(dar_folder, paste0(contrast_name, ".csv")))
}

# Save and annotate differential accessibility
dar.df <- reduce(dar.l, full_join, by = "peaks")
write_csv(dar.df, file.path(dar_folder, "all_comparisons.csv"))

for (name in names(dar.l)) {
  df <- dar.l[[name]]
  logFC_col <- paste0("logFC_", name)
  FDR_col <- paste0("FDR_", name)
  
  df <- df %>%
    mutate(diffaccessible = case_when(
      !!sym(FDR_col) <= 0.05 & !!sym(logFC_col) >= 1 ~ "Up",
      !!sym(FDR_col) <= 0.05 & !!sym(logFC_col) <= -1 ~ "Down",
      TRUE ~ "No")) %>%
    arrange(desc(!!sym(logFC_col)))
  
  write_csv(df, file.path(dar_difexp_folder, paste0(name, "_diffacc_full.csv")))
  write_csv(filter(df, diffaccessible != "No"), file.path(dar_difexp_folder, paste0(name, "_diffacc.csv")))
}


# =========================
# Volcano plots (all contrasts)
# =========================


volcano_dir <- "../results/plots/graphs_ATAC"
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

# Display name mapping
full_names <- c("Artery" = "Artery",
                "CapArt" = "Capillary artery",
                "Cap" = "Capillary",
                "CapVein" = "Capillary vein",
                "LargeVein" = "Large vein")

# Loop and plot
for (name in names(dar.l)) {
  df <- dar.l[[name]]
  
  logFC_col <- paste0("logFC_", name)
  FDR_col <- paste0("FDR_", name)
  
  df <- df %>%
    column_to_rownames("peaks") %>%
    arrange(desc(!!sym(logFC_col))) %>%
    mutate(diffaccesible = case_when(!!sym(FDR_col) <= 0.05 & !!sym(logFC_col) >= 1 ~ "Up",
                                     !!sym(FDR_col) <= 0.05 & !!sym(logFC_col) <= -1 ~ "Down",
                                     TRUE ~ "No"))
  
  summary_data <- df %>%
    group_by(diffaccesible) %>%
    summarise(n = n(), .groups = "drop")
  
  abbrev <- stringr::str_extract(name, "^[A-Za-z]+")
  full_name <- full_names[abbrev]
  subtitle <- glue::glue("{summary_data$n[summary_data$diffaccesible == 'Down']} closed regions and {summary_data$n[summary_data$diffaccesible == 'Up']} opened regions")
  comparison_type <- ifelse(stringr::str_detect(name, "KO_vs_Flox"), 
                            "<br><b>DARs </b><b><i>Pdcd10</b></i><sup><b>BECKO</b></sup> vs <b><i>Pdcd10</b></i><sup><b>fl/fl</b></sup><b></b>",
                            "<br><b>DARs vs </b><i><b>Pdcd10</b></i><sup><b>fl/fl</b></sup><b> all BEC subtypes</b>")
  
  volcano_p <- EnhancedVolcano(df,
                               lab = NA,
                               x = logFC_col,
                               y = FDR_col,
                               title = glue::glue("{full_name} {comparison_type}"),
                               subtitle = subtitle,
                               subtitleLabSize = 15,
                               pCutoff = 0.05,
                               FCcutoff = 1,
                               legendLabSize = 12,
                               axisLabSize = 13) + theme(plot.title = element_markdown())
  
  ggsave(filename = file.path(volcano_dir, paste0(name, "_volcano.png")),
         plot = volcano_p, width = 6, height = 5, dpi = 300)
}

# =========================
# DEG summary barplots
# =========================


summary_dir <- "../results/plots/graphs_endo"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# Load differential expression files
dar_files <- list.files("../results/ATAC/edgeR/dar_diffacc_tables/", pattern = "_diffacc.csv", full.names = TRUE)
samples <- basename(dar_files) %>% str_remove("_diffacc.csv")

# Build summary table
peaks_df_l <- map2(dar_files, samples, ~ {
  read.csv(.x) %>%
    select(peaks, diffaccessible) %>%
    mutate(cell_type = str_replace(.y, "_vs.*", ""))
})

peaks_df <- bind_rows(peaks_df_l)

peak_summ <- peaks_df %>%
  count(cell_type, diffaccessible, name = "num_peaks") %>%
  mutate(condition = if_else(str_detect(cell_type, "KO"), "KO", "Flox"),
         cell_type = str_replace_all(cell_type, 
                                     c("_KO" = "",
                                       "CapArt" = "Capillary artery", 
                                       "CapVein|CapVen" = "Capillary vein",
                                       "LargeVein" = "Large vein", 
                                       "\\bCap\\b" = "Capillary")),
         diffaccessible = if_else(diffaccessible == "Up", "Accesible\ncCRE", "Non-accesible\ncCRE"))

write.csv(peak_summ, "../results/ATAC/data/summary_DARs.csv", row.names = FALSE)

# Color palette
EC_colors <- rev(c("#247567", "#E2705B", "#6C7A99", "#B4698B", "#D8A21F"))

# -------- FLOX
peak_summ_Fl <- filter(peak_summ, condition == "Flox")
peak_summ_Fl$cell_type <- factor(peak_summ_Fl$cell_type, levels = rev(c("Artery", "Capillary artery",
                                                                        "Capillary", "Capillary vein",
                                                                        "Large vein")))

peak_fl_p <- ggplot(peak_summ_Fl, aes(y = cell_type)) +
  geom_segment(aes(x = 0,
                   xend = ifelse(diffaccessible == "Non-accesible\ncCRE", -num_peaks, num_peaks),
                   yend = cell_type,
                   color = diffaccessible), linewidth = 1) +
  geom_point(aes(x = ifelse(diffaccessible == "Non-accesible\ncCRE", -num_peaks, num_peaks),
                 color = diffaccessible), size = 3) +
  geom_text(aes(x = ifelse(diffaccessible == "Non-accesible\ncCRE", -num_peaks, num_peaks),
                label = abs(num_peaks),
                hjust = ifelse(diffaccessible == "Non-accesible\ncCRE", 1.23, -0.3)),
            size = 3.5) +
  scale_color_manual(values = c("Non-accesible\ncCRE" = "dodgerblue", "Accesible\ncCRE" = "#EE6F5D")) +
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_text(size = 11, margin = margin(t = 8, unit = "pt")),
        axis.text.y = element_text(size = 10, hjust = 1, face = "bold", colour = EC_colors),
        legend.text = element_text(size = 9.5, face = "bold", margin = margin(b = 5)),
        axis.title.y = element_blank(),
        plot.title = ggtext::element_markdown(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9)) +
  labs(x = "Number of significant DARs", y = "", color = "",
       title = "<b>BEC subtype marker cCRE</b>
         <b>(</b><i><b>Pdcd10</b></i><sup><b>fl/fl</b></sup><b>)</b>",
       subtitle = "Accesible and non-accesible chromatin |FC| > 1") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25) +
  scale_x_continuous(breaks = seq(-15000, 15000, by = 5000), labels = function(x) abs(x), limits = c(-13000, 16000)) 


ggsave(file.path(summary_dir, "lollipop_peaks_Fl.png"), peak_fl_p, width = 5, height = 3, dpi = 500)

# -------- KO
peak_summ_KO <- filter(peak_summ, condition == "KO")
peak_summ_KO$cell_type <- paste(peak_summ_KO$cell_type, "KO")
peak_summ_KO$cell_type <- factor(peak_summ_KO$cell_type, levels = rev(paste(c("Artery", "Capillary artery",
                                                                              "Capillary", "Capillary vein",
                                                                              "Large vein"), "KO")))

peak_ko_p <- ggplot(peak_summ_KO, aes(y = cell_type)) +
  geom_segment(aes(x = 0,
                   xend = ifelse(diffaccessible == "Non-accesible\ncCRE", -num_peaks, num_peaks),
                   yend = cell_type,
                   color = diffaccessible), linewidth = 1) +
  geom_point(aes(x = ifelse(diffaccessible == "Non-accesible\ncCRE", -num_peaks, num_peaks),
                 color = diffaccessible), size = 3) +
  geom_text(aes(x = ifelse(diffaccessible == "Non-accesible\ncCRE", -num_peaks, num_peaks),
                label = abs(num_peaks),
                hjust = ifelse(diffaccessible == "Non-accesible\ncCRE", 1.23, -0.2)),
            size = 3.5) +
  scale_color_manual(values = c("Non-accesible\ncCRE" = "dodgerblue", "Accesible\ncCRE" = "#EE6F5D")) +
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_text(size = 11, margin = margin(t = 8, unit = "pt")),
        axis.text.y = element_text(size = 10, hjust = 1, face = "bold", colour = EC_colors),
        legend.text = element_text(size = 9, face = "bold", margin = margin(b = 5)),
        axis.title.y = element_blank(),
        plot.title = ggtext::element_markdown(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9)) +
  labs(x = "Number of significant DARs", y = "", color = "",
       title = "<b>DA cCRE per BEC subtype</b><br>
             <b>(</b><i><b>PDCD10</b></i><sup><b>BECKO</b></sup> 
             <b> vs </b> <i><b>PDCD10</b></i><sup><b>fl/fl</b></sup><b>)</b>",
       subtitle = "Accesible and non-accesible chromatin |FC| > 1") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25) +
  scale_x_continuous(breaks =  seq(-30000, 30000, by = 15000), labels = function(x) abs(x), limits = c(-18500, 32000)) 


ggsave(file.path(summary_dir, "lollipop_peaks_KO.png"), peak_ko_p, width = 4.5, height = 3, dpi = 500)

# =========================
# Upset plots KO
# =========================

# Load differential expression files
dar_files <- list.files("../results/ATAC/edgeR/dar_diffacc_tables/", pattern = "_vs_Flox_diffacc.csv", full.names = TRUE)
samples <- basename(dar_files) %>% str_remove("_vs_Flox_diffacc.csv")

# Build summary table
peaks_df_l <- map2(dar_files, samples, ~ {
  read.csv(.x) %>%
    filter(diffaccessible == "Up") %>% pull("peaks") %>% unique()
})

peaks_df_l <- peaks_df_l[c(1, 3, 2, 4, 5)]

names(peaks_df_l) <- c("Artery", "Capillary artery","Capillary","Capillary vein","Large vein")

from_list <- function(list_data) {
  members = unique(unlist(list_data))
  data.frame(
    lapply(list_data, function(set) members %in% set),
    row.names=members,
    check.names=FALSE
  )
}


set_order <- rev(c("Artery", "Capillary artery","Capillary","Capillary vein","Large vein"))

peak_up_mtx <- from_list(peaks_df_l)

upset_peak_up <- upset(peak_up_mtx, set_order, n_intersections = 10, 
                       keep_empty_groups = T, sort_sets = F, name = NULL,
                       base_annotations = list('Intersection size'= (intersection_size(text = list(size = 4.5, vjust = -0.2),
                                                                                       bar_number_threshold = 1,
                                                                                       fill='grey80') ) +
                                                 labs(title = "<b>Intersection of gained DARs</b><br>
                                                <b>(</b><i><b>PDCD10</b></i><sup><b>BECKO</b></sup> 
                                                     <b> vs </b> <i><b>PDCD10</b></i><sup><b>fl/fl</b></sup><b>)</b>") + 
                                                 expand_limits(y = 1300) +
                                                 scale_y_continuous(limits = c(0, 5250)) +
                                                 theme(text = element_text(size = 12), 
                                                       axis.text.y = element_text(size = 14),
                                                       axis.title.y = element_text(size = 15, vjust = -20, margin = margin(r = 20)),
                                                       plot.title = ggtext::element_markdown(hjust = 0.5, size = 12))),
                       height_ratio = 0.3, set_sizes = F) + theme(axis.text.y = element_text(size = 12, face = "bold")) 


ggsave(filename = "../results/plots/graphs_endo/upset_peaks_up_KO.png", plot = upset_peak_up, width = 7, height = 5, bg = "white")

# Optionally save session info
sink("../results/logs/sessioninfo_09_edgeR_diff_accesibility.txt")
sessionInfo()
sink()
