pacman::p_load( "limma", "tibble", "viridis", "tidyr", "readr",
                "RColorBrewer", "pheatmap", "circlize", "stringr",
                "dplyr", "ggplot2", "textshape", "edgeR", "paletteer")

### ATAC ----

# Create a folder to save the CSV files

dar_folder <- "../results/Human_cells/D3_ATAC/edgeR/dar_tables"
dar_difacc_tables <- "../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables"
graphs <- "../results/Human_cells/D3_ATAC/edgeR/plots"

# graph volcano folder

dir.create(dar_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(dar_difacc_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(graphs, recursive = TRUE, showWarnings = FALSE)

comparison_names <- c("KD_VEH_vs_WT_VEH", "KD_VEH_vs_WT_CCM",
                      "KD_CCM_vs_KD_VEH", "KD_CCM_vs_WT_VEH", "KD_CCM_vs_WT_CCM",
                      "KD_T52_vs_KD_VEH", "KD_T52_vs_WT_VEH", "KD_T52_vs_WT_T52",
                      "KD_CCT_vs_KD_CCM", "KD_CCT_vs_KD_T52", "KD_CCT_vs_KD_VEH",
                      "KD_CCT_vs_WT_VEH", "KD_CCT_vs_WT_T52", "KD_CCT_vs_WT_CCT",
                      "WT_CCM_vs_WT_VEH", "WT_T52_vs_WT_VEH", "WT_CCT_vs_WT_VEH",
                      "K1_CCM_vs_K1_VEH", "K1_CCM_vs_W1_VEH", "W1_CCM_vs_W1_VEH", 
                      "KD_CCM_vs_K1_CCM", "WT_CCM_vs_W1_CCM")


###union peaksets
union_peakset <- read.table("../results/Human_cells/D3_ATAC/D3_merged_peaks/celltype_condition.filteredNfixed.union.peakSet", sep = "\t", header = TRUE) %>% 
  filter(spm > 2) %>% 
  mutate("PeakID" = paste(seqnames, start, end, sep = "-"))

peaks_tofilt <- union_peakset$PeakID

files <- list.files("../results/Human_cells/D3_ATAC/D3_merged_peaks/peakset_by_comparisons/")


samples <- gsub(".filteredNfixed.union.peakset", "", files)[c(17, 16, 10, 12, 11, 13, 15, 14, 3, 4, 5, 8, 7, 6, 21, 22, 19, 1, 2, 18, 9, 20)]

names(comparison_names) <- samples

peaks_by_comp <- list()
for (sample in samples) {
  
  comparison_name <- comparison_names[[sample]]
  
  tmp <-  read.table(sprintf("../results/Human_cells/D3_ATAC/D3_merged_peaks/peakset_by_comparisons/%s.filteredNfixed.union.peakset", sample), sep = "\t", header = TRUE) %>% 
    mutate("PeakID" = paste(chr, start, end, sep = "-")) %>%  
    filter(PeakID %in% peaks_tofilt) %>% 
    pull("PeakID")
  
  peaks_by_comp[[comparison_name]] <- tmp
  
  rm(tmp, sample, comparison_name)
}

##Counts----

cts <- read.csv("../results/Human_cells/D3_ATAC/data/cts_batch_corrected.csv", row.names = 1)

cpm_atac <- cpm(cts)
cpm_atac <- as.data.frame(cpm_atac)

log2cpm_atac <- log2(cpm_atac + 1)
log2cpm_atac <- as.data.frame(log2cpm_atac)

cpm_atac <- rownames_to_column(.data = cpm_atac, var = "peaks")
write_csv(x = cpm_atac, file = "../results/Human_cells/D3_ATAC/data/cpm.csv")


log2cpm_atac <- rownames_to_column(.data = log2cpm_atac, var = "peaks")
write.csv(x = log2cpm_atac, file = "../results/Human_cells/D3_ATAC/data/log2cpm.csv")

#groups

# Define the conditions and their replicates (3 replicates each)
replicates <- 3

# Generate the list based on the comparisons
groups_cts_l <- list()
for (comp in comparison_names) {
  conditions <- strsplit(comp, "_vs_")[[1]]
  groups_cts_l[[comp]] <- rep(conditions, each = replicates)
}

contrasts_v <- c("KD_VEH_vs_WT_VEH" = "KD_VEH - WT_VEH",
                 "KD_VEH_vs_WT_CCM" = "KD_VEH - WT_CCM",
                 "KD_CCM_vs_KD_VEH" = "KD_CCM - KD_VEH",
                 "KD_CCM_vs_WT_VEH" = "KD_CCM - WT_VEH",
                 "KD_CCM_vs_WT_CCM" = "KD_CCM - WT_CCM",
                 "KD_T52_vs_KD_VEH" = "KD_T52 - KD_VEH",
                 "KD_T52_vs_WT_VEH" = "KD_T52 - WT_VEH",
                 "KD_T52_vs_WT_T52" = "KD_T52 - WT_T52",
                 "KD_CCT_vs_KD_CCM" = "KD_CCT - KD_CCM",
                 "KD_CCT_vs_KD_T52" = "KD_CCT - KD_T52",
                 "KD_CCT_vs_KD_VEH" = "KD_CCT - KD_VEH",
                 "KD_CCT_vs_WT_VEH" = "KD_CCT - WT_VEH",
                 "KD_CCT_vs_WT_T52" = "KD_CCT - WT_T52",
                 "KD_CCT_vs_WT_CCT" = "KD_CCT - WT_CCT",
                 "WT_CCM_vs_WT_VEH" = "WT_CCM - WT_VEH",
                 "WT_T52_vs_WT_VEH" = "WT_T52 - WT_VEH",
                 "WT_CCT_vs_WT_VEH" = "WT_CCT - WT_VEH",
                 "K1_CCM_vs_K1_VEH" = "K1_CCM - K1_VEH", 
                 "K1_CCM_vs_W1_VEH" = "K1_CCM - W1_VEH",
                 "W1_CCM_vs_W1_VEH" = "W1_CCM - W1_VEH",
                 "KD_CCM_vs_K1_CCM" = "KD_CCM - K1_CCM",
                 "WT_CCM_vs_W1_CCM" = "WT_CCM - W1_CCM")

contrast_numbers <- list(KD_VEH_vs_WT_VEH = 1,
                         KD_VEH_vs_WT_CCM = 2, 
                         KD_CCM_vs_KD_VEH = 3,
                         KD_CCM_vs_WT_VEH = 4,
                         KD_CCM_vs_WT_CCM = 5,
                         KD_T52_vs_KD_VEH = 6,
                         KD_T52_vs_WT_VEH = 7,
                         KD_T52_vs_WT_T52 = 8,
                         KD_CCT_vs_KD_CCM = 9,
                         KD_CCT_vs_KD_T52 = 10,
                         KD_CCT_vs_KD_VEH = 11,
                         KD_CCT_vs_WT_VEH = 12,
                         KD_CCT_vs_WT_T52 = 13,
                         KD_CCT_vs_WT_CCT = 14,
                         WT_CCM_vs_WT_VEH = 15,
                         WT_T52_vs_WT_VEH = 16,
                         WT_CCT_vs_WT_VEH = 17,
                         K1_CCM_vs_K1_VEH = 18,
                         K1_CCM_vs_W1_VEH = 19,
                         W1_CCM_vs_W1_VEH = 20, 
                         KD_CCM_vs_K1_CCM = 21,
                         WT_CCM_vs_W1_CCM = 22)

cts_matches <- c("KD_V1|KD_V2|KD_V3|WT_V1|WT_V2|WT_V3",  # KD_VEH_vs_WT_VEH
                 "KD_V1|KD_V2|KD_V3|WT_C1|WT_C2|WT_C3",  # KD_VEH_vs_WT_CCM
                 "KD_C1|KD_C2|KD_C3|KD_V1|KD_V2|KD_V3",  # KD_CCM_vs_KD_VEH
                 "KD_C1|KD_C2|KD_C3|WT_V1|WT_V2|WT_V3",  # KD_CCM_vs_WT_VEH
                 "KD_C1|KD_C2|KD_C3|WT_C1|WT_C2|WT_C3",  # KD_CCM_vs_WT_CCM
                 "KD_T1|KD_T2|KD_T3|KD_V1|KD_V2|KD_V3",  # KD_T52_vs_KD_VEH
                 "KD_T1|KD_T2|KD_T3|WT_V1|WT_V2|WT_V3",  # KD_T52_vs_WT_VEH
                 "KD_T1|KD_T2|KD_T3|WT_T1|WT_T2|WT_T3",  # KD_T52_vs_WT_T52
                 "KC_T1|KC_T2|KC_T3|KD_C1|KD_C2|KD_C3",  # KD_CCT_vs_KD_CCM
                 "KC_T1|KC_T2|KC_T3|KD_T1|KD_T2|KD_T3",  # KD_CCT_vs_KD_T52
                 "KC_T1|KC_T2|KC_T3|KD_V1|KD_V2|KD_V3",  # KD_CCT_vs_KD_VEH
                 "KC_T1|KC_T2|KC_T3|WT_V1|WT_V2|WT_V3",  # KD_CCT_vs_WT_VEH
                 "KC_T1|KC_T2|KC_T3|WT_T1|WT_T2|WT_T3",  # KD_CCT_vs_WT_T52
                 "KC_T1|KC_T2|KC_T3|WC_T1|WC_T2|WC_T3",  # KD_CCT_vs_WT_CCT
                 "WT_C1|WT_C2|WT_C3|WT_V1|WT_V2|WT_V3",  # WT_CCM_vs_WT_VEH
                 "WT_T1|WT_T2|WT_T3|WT_V1|WT_V2|WT_V3",  # WT_T52_vs_WT_VEH
                 "WC_T1|WC_T2|WC_T3|WT_V1|WT_V2|WT_V3",  # WT_CCT_vs_WT_VEH 
                 "K1_C1|K1_C2|K1_C3|K1_V1|K1_V2|K1_V3",  # K1_CCM_vs_K1_VEH
                 "K1_C1|K1_C2|K1_C3|W1_V1|W1_V2|W1_V3",  # K1_CCM_vs_W1_VEH
                 "W1_C1|W1_C2|W1_C3|W1_V1|W1_V2|W1_V3",  # W1_CCM_vs_W1_VEH 
                 "KD_C1|KD_C2|KD_C3|K1_C1|K1_C2|K1_C3",  # KD_CCM_vs_K1_CCM
                 "WT_C1|WT_C2|WT_C3|W1_C1|W1_C2|W1_C3")  # WT_CCM_vs_W1_CCM)

names(cts_matches) <- names(groups_cts_l)


peak.dge.celltype.l <- list()
peak.contrasts.l <- list()
for (loop_var in names(groups_cts_l)) {
  
  print(paste0("working on ", loop_var))
  
  cts_match <- cts_matches[[loop_var]]
  
  cts_tmp <- cts %>% 
    select(matches(cts_match)) 
 
  cts_tmp <- cts_tmp[, order(match(colnames(cts_tmp), unlist(strsplit(cts_match, "\\|"))))]
  
  print(colnames(cts_tmp))
  
  peaks_comp <- peaks_by_comp[[loop_var]]
  
  print(paste0("working on cpm ", loop_var))
  
  cpm_atac_tmp <- cpm(cts_tmp)
  cpm_atac_tmp <- as.data.frame(cpm_atac_tmp)

  cpm_atac_filt <- rownames_to_column(.data = cpm_atac_tmp, var = "peaks") %>% 
    filter(peaks %in% peaks_comp)
  
  print(paste0("working on selecting groups for ", loop_var))
  ctsGroups <- groups_cts_l[[loop_var]]
  
  
  d <- DGEList(counts = cts_tmp,
               group = factor(ctsGroups))
  
  d
  
  tmp_peaks <- peaks_comp[peaks_comp %in% rownames(d)]
  
  
  d <- d[tmp_peaks, ]
  
  d$samples$lib.size <- colSums(d$counts)
  
  d <- calcNormFactors(d)
  
  png(filename = paste0(graphs, "/", loop_var, "_MDSplot.png"), width = 6, height = 4, units = "in", res = 300)
  #plotMDS(d, method="bcv", col=as.numeric(d$samples$group), pch=20)
  plotMDS(d,col=as.numeric(d$samples$group), pch=20)
  legend("center", as.character(unique(d$samples$group)), col=1:6, pch=20)
  dev.off()

  design.mat <- model.matrix(~ 0 + d$samples$group)
  colnames(design.mat) <- levels(d$samples$group)
  
  d <- estimateDisp(d, design.mat)
  
  # Fit the negative binomial GLM
  fit <- glmFit(d, design.mat)
  
  print(paste0("working on selecting contrasts ", loop_var))
  contrast_tmp <- contrast_numbers[[loop_var]]
  
  contrasts <- makeContrasts(contrasts_v[contrast_tmp], levels = design.mat)
  colnames(contrasts) <- names(contrasts_v[contrast_tmp])
  
  
  ### Fit the contrasts loop for each cell type
  
  # Number of contrasts
  
  contrast_name <- colnames(contrasts)
  
  # Perform LRT for the current contrast
  fit_contrast <- glmLRT(fit, contrast = contrasts)
  
  
  peak.contrasts.l[[contrast_name]] <- fit_contrast
  peak.contrasts.l[[contrast_name]]$comparison <- paste0(contrast_name)
  
  print(paste0("working on extracting results ", loop_var))
  # Extract top tags
  results_table <- topTags(fit_contrast, adjust.method = "fdr", n = Inf)$table 
  results_table <- results_table %>% 
    rename_with(~ paste0(., "_", contrast_name)) %>% rownames_to_column("peaks") %>% 
    left_join(cpm_atac_filt, by = "peaks")
  
  
  peak.dge.celltype.l[[contrast_name]] <- results_table
  
  # Save results to a CSV file
  file_name <- paste0(dar_folder, "/", contrast_name, ".csv")
  write_csv(results_table, file = file_name)

}


dge_values <- names(peak.dge.celltype.l)

for (loop_var in dge_values) {
  data <- peak.dge.celltype.l[[loop_var]]
  data_diff <- data %>% 
    mutate(diffaccessible = ifelse(data[ ,6] <= 0.05 & data[ ,2] >= 1, "Up", 
                                  ifelse(data[ ,6] <= 0.05 & data[ ,2] <= -1, "Down", "No"))) %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var))))
  
  
  data_diff_filt <- data_diff %>% 
    dplyr::filter(!diffaccessible == "No")
  
  write.csv(x = data_diff_filt, file = sprintf("../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables/%s_diffacc.csv", loop_var), row.names = F)
  write.csv(x = data_diff, file = sprintf("../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables/%s_diffacc_full.csv", loop_var), row.names = F)
  
}


files <- list.files(path = "../results/Human_cells/D3_ATAC/edgeR/dar_tables/", recursive = F)
files 

samples <- gsub(".csv", "", files)
samples

#volcano plots----
data_volcano <- list()
for (loop_var in samples) {
  data <- read.csv(file = sprintf("../results/Human_cells/D3_ATAC/edgeR/dar_tables/%s.csv",loop_var))
  data <- data %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var)))) %>% 
    mutate(diffaccessible = ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) >= 1, "Up", 
                                  ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) <= -1, "Down", "No"))) %>% 
    dplyr::select(peaks, starts_with("logFC"), starts_with("FDR"), diffaccessible) %>% 
    column_to_rownames("peaks")
  
  data_volcano[[loop_var]] <- data
  
}


full_names <- str_replace_all(samples, rev(c("KD" = "siPDCD10","WT" = "Control", "K1" = "siPDCD10 1 hour", "W1" = "Control 1 hour", "CCT" = "(CCM-like env + T5224)", 
                                             "T52" = "(T5224)", "CCM" = "(CCM-like env)", 
                                             "VEH" = "(Vehicle)", "_" = " ")))

names(full_names) <- samples
full_names

volcano_plots <- list()
for (loop_var in names(data_volcano)) {
  data <- data_volcano[[loop_var]]
  summary_data <- data %>% 
    group_by(diffaccessible) %>% 
    summarise(n = n())
  
  full_name <- full_names[loop_var]
  
  title <- paste0("DAR\n", full_name)
  
  volcano_p <- EnhancedVolcano::EnhancedVolcano(data, 
                                                lab = NA,
                                                titleLabSize = 14,
                                                subtitleLabSize = 12,
                                                axisLabSize = 12,
                                                legendLabSize = 11,
                                                legendIconSize = 3.5,
                                                x = sprintf("logFC_%s", loop_var),
                                                y = sprintf("FDR_%s", loop_var), 
                                                title = title, 
                                                subtitle = paste0(summary_data[1,]$n, " accesible DARs and ", summary_data[3,]$n, " non-accesible DARs"),
                                                pCutoff = 0.05,
                                                FCcutoff = 1, 
                                                raster = T) 
  ggsave(filename = paste0(graphs, "/", loop_var, "_volcano.png"),
         plot = volcano_p, width = 6, height = 6, units = "in", dpi = 300, device = "png")
  
}


# Optionally save session info
sink("../logs/sessioninfo_human_10_DAR_edgeR.txt")
sessionInfo()
sink()









