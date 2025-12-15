
pacman::p_load( "limma", "tibble", "viridis", "tidyr",
                "RColorBrewer", "pheatmap", "circlize", "stringr",
                "dplyr", "ggplot2", "textshape", "edgeR", "paletteer")

### RNA ----

# Create a folder to save the CSV files

deg_folder <- "../results/Human_cells/D3_RNA/edgeR/deg_tables"
if (!dir.exists(deg_folder)){
  dir.create(deg_folder, recursive = TRUE)
}

deg_difexp_folder <- "../results/Human_cells/D3_RNA/edgeR/deg_difexp_tables"
if (!dir.exists(deg_difexp_folder)){
  dir.create(deg_difexp_folder, recursive = TRUE)
}

data_folder <- "../results/Human_cells/D3_RNA/data"
if (!dir.exists(data_folder)){
  dir.create(data_folder, recursive = TRUE)
} 


# graph folder
if (!dir.exists("../results/Human_cells/D3_RNA/graphs_RNA")){
  dir.create("../results/Human_cells/D3_RNA/graphs_RNA")
}

# graph MDS folder
if (!dir.exists("../results/Human_cells/D3_RNA/graphs_RNA/MDS")){
  dir.create("../results/Human_cells/D3_RNA/graphs_RNA/MDS", recursive = T)
}

# graph volcano folder
if (!dir.exists("../results/Human_cells/D3_RNA/graphs_RNA/volcano")){
  dir.create("../results/Human_cells/D3_RNA/graphs_RNA/volcano", recursive = T)
}

# Load raw count data

cts <- read.table(file = "../results/Human_cells/D3_RNA/data/salmon.merged.gene_counts_length_scaled.tsv", sep = "\t", header = TRUE) %>%  
  dplyr::select(-gene_id)

colnames(cts)

cts <- cts %>% 
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames("gene_name") %>% 
  dplyr::select(starts_with("WT_V"), starts_with("KD_V"),
         starts_with("WT_C"), starts_with("KD_C"),
         starts_with("WT_T"), starts_with("KD_T"),
         starts_with("WC_T"), starts_with("KC_T"))

#write.csv(cts, file = "data/counts_ordered.csv")

exprdata <- as.matrix(cts)

metadf <- data.frame(Batch = c(rep(c(rep(c("A", "B"), c(2, 1))), 4), rep("B", 12) ), 
                     Group = rep(c(rep(c("WT", "KD"), c(3,3))), 4),
                     treatment = rep(c("veh", "ccm",  "t52", "ct52"), c(6, 6,  6, 6)),
                     row.names = colnames(exprdata))


table(metadf$Batch, interaction(metadf$Group, metadf$treatment, drop = TRUE))

# cols for AB cohort (veh/ccm only)
ab_cols <- rownames(metadf)[metadf$Batch %in% c("A","B") &
                              metadf$treatment %in% c("veh","ccm")]

cts_ab   <- as.matrix(cts[, ab_cols])
meta_ab  <- metadf[ab_cols, , drop = FALSE]

library(sva)
cb_counts_ab <- ComBat_seq(counts = cts_ab,
                           batch  = meta_ab$Batch,
                           group  = interaction(meta_ab$Group, meta_ab$treatment, drop = TRUE))


cts <- as.data.frame(cbind(cb_counts_ab, exprdata[ ,13:24]))

cpm_RNA <- cpm(cts)
cpm_RNA <- as.data.frame(cpm_RNA)

log2cpm_RNA <- log2(cpm_RNA + 1)
log2cpm_RNA <- as.data.frame(log2cpm_RNA)

cpm_RNA <- rownames_to_column(.data = cpm_RNA, var = "gene")
#write.csv(x = cpm_RNA, file = "../results/Human_cells/D3_RNA/data/cpm.csv", row.names = F)


log2cpm_RNA <- rownames_to_column(.data = log2cpm_RNA, var = "gene")
# write.csv(x = log2cpm_RNA, file = "../results/Human_cells/D3_RNA/data/log2cpm.csv", row.names = F)


avg_logCPM <- rowMeans(log2cpm_RNA[ ,-1])

keep <- avg_logCPM > 1.5
cpm_RNA_f <- cpm_RNA[keep,]
log2cpm_RNA_f <- log2cpm_RNA[keep,]
write.csv(x = cpm_RNA_f, file = "../results/Human_cells/D3_RNA/data/cpm_filt.csv", row.names = F)
write.csv(x = log2cpm_RNA_f, file = "../results/Human_cells/D3_RNA/data/log2cpm_filt.csv", row.names = F)

##
comparison_names <- c("KD_VEH_vs_WT_VEH", "KD_VEH_vs_WT_CCM",
                      "KD_CCM_vs_KD_VEH", "KD_CCM_vs_WT_VEH", "KD_CCM_vs_WT_CCM",
                      "KD_T52_vs_KD_VEH", "KD_T52_vs_WT_VEH", "KD_T52_vs_WT_T52",
                      "KD_CCT_vs_KD_CCM", "KD_CCT_vs_KD_T52", "KD_CCT_vs_KD_VEH",
                      "KD_CCT_vs_WT_VEH", "KD_CCT_vs_WT_T52", "KD_CCT_vs_WT_CCT",
                      "WT_CCM_vs_WT_VEH", "WT_T52_vs_WT_VEH", "WT_CCT_vs_WT_VEH")

# Define the conditions and their replicates (3 replicates each)
replicates <- 3

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
                 "WT_CCT_vs_WT_VEH" = "WT_CCT - WT_VEH")

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
                         WT_CCT_vs_WT_VEH = 17)

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
                 "WC_T1|WC_T2|WC_T3|WT_V1|WT_V2|WT_V3")   # WT_CCT_vs_WT_VEH 

names(cts_matches) <- names(groups_cts_l)



gene.dge.celltype.l <- list()
gene.contrasts.l <- list()
for (loop_var in names(groups_cts_l)) {
  
  print(paste0("working on ", loop_var))
  
  cts_match <- cts_matches[[loop_var]]
  
  cts_tmp <- cts %>% 
    dplyr::select(matches(cts_match)) 
  
  cts_tmp <- cts_tmp[, order(match(colnames(cts_tmp), unlist(strsplit(cts_match, "\\|"))))]
  
  print(colnames(cts_tmp))
  
  
  print(paste0("working on cpm ", loop_var))
  
  cpm_atac <- cpm(cts_tmp)
  cpm_atac <- as.data.frame(cpm_atac)
  
  log2cpm_atac <- log2(cpm_atac + 1)
  log2cpm_atac <- as.data.frame(log2cpm_atac)
  
  cpm_atac <- rownames_to_column(.data = cpm_atac, var = "gene") 
  
  #saveRDS(object = cpm_atac, file = paste0("rds_bulk_RNA/cpm_ko_", loop_var, ".rds"))
  
  log2cpm_atac <- rownames_to_column(.data = log2cpm_atac, var = "gene") 
  
  #saveRDS(object = log2cpm_atac, file = paste0("rds_bulk_RNA/log2cpm_ko_", loop_var, ".rds"))
  
  print(paste0("working on selecting groups for ", loop_var))
  ctsGroups <- groups_cts_l[[loop_var]]
  
  
  d <- DGEList(counts = cts_tmp,
               group = factor(ctsGroups))
  
  d

  keep <- avg_logCPM > 1.5  # Adjust this based on your histogram
  d <- d[keep,]
  
  d$samples$lib.size <- colSums(d$counts)
  
  d <- calcNormFactors(d)
  
  
  png(filename = paste0("../results/Human_cells/D3_RNA/graphs_RNA/MDS/", loop_var, "_MDSplot.png"), width = 6, height = 4, units = "in", res = 300)
  # Adjust margins for space on the right for the legend
  par(mar = c(5, 4, 4, 8))  # Increase the right margin (4th value)
  # Plot MDS
  plotMDS(d, col = as.numeric(d$samples$group), pch = 20)
  # Add the legend outside the plot area
  legend("topright", inset = c(-0.3, 0), legend = as.character(unique(d$samples$group)), col = 1:6, pch = 20, xpd = TRUE)
  dev.off()
  
  
  design.mat <- model.matrix(~ 0 + d$samples$group)
  colnames(design.mat) <- levels(d$samples$group)
  
  d1 <- estimateDisp(d, design.mat)
  dim(d1$counts)
  
  # Fit the negative binomial GLM
  fit <- glmFit(d1, design.mat)
  
  print(paste0("working on selecting contrasts ", loop_var))
  contrast_tmp <- contrast_numbers[[loop_var]]
  
  contrasts <- makeContrasts(contrasts_v[contrast_tmp], levels = design.mat)
  colnames(contrasts) <- names(contrasts_v[contrast_tmp])
  
  
  ### Fit the contrasts loop for each cell type
  
  # Number of contrasts
  
  contrast_name <- colnames(contrasts)
  
  # Perform LRT for the current contrast
  fit_contrast <- glmLRT(fit, contrast = contrasts)
  
  
  gene.contrasts.l[[contrast_name]] <- fit_contrast
  gene.contrasts.l[[contrast_name]]$comparison <- paste0(contrast_name)
  
  print(paste0("working on extracting results ", loop_var))
  # Extract top tags
  results_table <- topTags(fit_contrast, adjust.method = "fdr", n = nrow(d1))$table 
  results_table <- results_table %>% 
    rename_all(~ paste0(., "_", contrast_name)) %>% rownames_to_column("gene") %>% 
    left_join(log2cpm_atac, by = "gene")
  
  
  gene.dge.celltype.l[[contrast_name]] <- results_table
  
  # Save results to a CSV file
  file_name <- paste0(deg_folder, "/", contrast_name, ".csv")
  write.csv(results_table, file = file_name, row.names = FALSE)
  # 
}

gene.dge.celltype.df <- purrr::reduce(gene.dge.celltype.l, full_join, by = "gene")

write.csv(gene.dge.celltype.df, 
          file = "../results/Human_cells/D3_RNA/edgeR/deg_tables/all_contrasts_RNA.csv",
          row.names = FALSE)


#The numbers of DE genes under each comparison are shown below.
dt <- lapply(lapply(gene.contrasts.l, decideTests), summary)
dt.all <- do.call("cbind", dt)
dt.all


dge_values <- names(gene.dge.celltype.l)

for (loop_var in dge_values) {
  data <- gene.dge.celltype.l[[loop_var]]
  data_diff <- data %>% 
    mutate(diffexpressed = ifelse(data[ ,6] <= 0.05 & data[ ,2] >= 0.5, "Up", 
                                  ifelse(data[ ,6] <= 0.05 & data[ ,2] <= -0.5, "Down", "No"))) %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var))))
  
  data_diff_filt <- data_diff %>% 
    dplyr::filter(!diffexpressed == "No")
  
  write.csv(x = data_diff_filt, file = sprintf("../results/Human_cells/D3_RNA/edgeR/deg_difexp_tables/%s_diffexp.csv", loop_var), row.names = F)
  write.csv(x = data_diff, file = sprintf("../results/Human_cells/D3_RNA/edgeR/deg_difexp_tables/%s_diffexp_full.csv", loop_var), row.names = F)
}


### volcanos

files <- list.files(path = "../results/Human_cells/D3_RNA/edgeR/deg_tables/", recursive = F)
files 

samples <- gsub(".csv", "", files[-1])


#volcano plots----

data_volcano <- list()
for (loop_var in samples) {
  data <- read.csv(file = sprintf("../results/Human_cells/D3_RNA/edgeR/deg_tables/%s.csv",loop_var))
  data <- data %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var)))) %>% 
    mutate(diffexpressed = ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) >= 0.5, "Up", 
                                  ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) <= -0.5, "Down", "No"))) %>% 
    dplyr::select(gene, starts_with("logFC"), starts_with("FDR"), diffexpressed) %>% 
    column_to_rownames("gene")
  
  data_volcano[[loop_var]] <- data
}

full_names <- str_replace_all(samples, rev(c("KD" = "siPDCD10","WT" = "Control", "CCT" = "(CCM-like env + T5224)", 
                                             "T52" = "(T5224)", "CCM" = "(CCM-like env)", 
                                             "VEH" = "(Vehicle)", "_" = " ")))

names(full_names) <- samples
full_names

volcano_plots <- list()
for (loop_var in names(data_volcano)) {
  data <- data_volcano[[loop_var]]
  summary_data <- data %>% 
    group_by(diffexpressed) %>% 
    summarise(n = n())
  
  full_name <- full_names[loop_var]
  
  title <- paste0("DEG\n", full_name)
  
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
                                                subtitle = paste0(summary_data[1,]$n, " downregulated and ", summary_data[3,]$n, " upregulated"),
                                                pCutoff = 0.05,
                                                FCcutoff = 0.5, 
                                                raster = T) 
  ggsave(filename = paste0("../results/Human_cells/D3_RNA/graphs_RNA/volcano/", loop_var, "_volcano.png"),
         plot = volcano_p, width = 6, height = 6, units = "in", dpi = 300, device = "png")
  #assign(paste0(loop_var, "_hmap"), value = hm)
  volcano_plots[[loop_var]] <- volcano_p
  
}

# Optionally save session info
sink("../logs/sessioninfo_human_11_DEG_EdgeR_batch_corrected.txt")
sessionInfo()
sink()



