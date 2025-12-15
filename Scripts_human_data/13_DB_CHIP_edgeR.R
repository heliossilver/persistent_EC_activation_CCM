
pacman::p_load( "limma", "tibble", "viridis", "tidyr",
                "RColorBrewer", "pheatmap", "circlize", "stringr",
                "dplyr", "ggplot2", "textshape", "edgeR", "paletteer")

### ChIP_JUNB ----

# Create a folder to save the CSV files

deg_folder <- "../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/deg_tables"
if (!dir.exists(deg_folder)){
  dir.create(deg_folder, recursive = TRUE)
}

deg_difexp_folder <- "../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/deg_difexp_tables"
if (!dir.exists(deg_difexp_folder)){
  dir.create(deg_difexp_folder, recursive = TRUE)
}


# graph folder
if (!dir.exists("../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/graphs_ChIP_JUNB")){
  dir.create("../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/graphs_ChIP_JUNB")
}


# graph volcano folder
if (!dir.exists("../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/graphs_ChIP_JUNB/volcano")){
  dir.create("../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/graphs_ChIP_JUNB/volcano", recursive = T)
}


###union peaksets
union_peakset_junb <- read.table("../results/Human_cells/ChipSeq/data/JUNB_mergedpeaks/JUNB.filteredNfixed.union.peakSet", sep = "\t")

colnames(union_peakset_junb) <- c("seqnames",	"start", "end",	"width",	"strand",	"score",	"name",	"label",	"GC",	"N",	"spm")

union_peakset_junb <- union_peakset_junb %>% 
  filter(spm > 2) %>% 
  mutate("PeakID" = paste(seqnames, start, end, sep = "-"))

peaks_tofilt <- union_peakset_junb$PeakID


files <- list.files("../results/Human_cells/ChipSeq/data/JUNB_mergedpeaks/peakset_by_comparisons/")
samples <- gsub(".filteredNfixed.union.peakset", "", files[c(5, 4, 1, 3, 2, 6)])


comparison_names <- c("KD_VEH_vs_WT_VEH", "KD_VEH_vs_WT_CCM",
                      "KD_CCM_vs_KD_VEH", "KD_CCM_vs_WT_VEH", 
                      "KD_CCM_vs_WT_CCM", "WT_CCM_vs_WT_VEH")

names(comparison_names) <- samples

peaks_by_comp <- list()
for (sample in samples) {
  
  comparison_name <- comparison_names[[sample]]
  
  tmp <-  read.table(sprintf("../results/Human_cells/ChipSeq/data/JUNB_mergedpeaks/peakset_by_comparisons/%s.filteredNfixed.union.peakset", sample), sep = "\t", header = TRUE) %>% 
    mutate("PeakID" = paste(chr, start, end, sep = "-")) %>%  
    filter(PeakID %in% peaks_tofilt) %>% 
    pull("PeakID")
  
  peaks_by_comp[[comparison_name]] <- tmp
  
  rm(tmp, sample, comparison_name)
}


cts <- read.table(file = "../results/Human_cells/ChipSeq/data/multicov_JUNB.tsv", sep = "\t")

colnames(cts) <- c("seqnames",	"start", "end",	"width",	"strand",	"score",	"name",	"label",	"GC",	"N",	"spm", 
                   "KD_C1", "KD_C2", "KD_V1", "KD_V2", "WT_C1", "WT_C2", "WT_V1", "WT_V2")

cts <-  cts %>% 
  mutate("PeakID" = paste(seqnames, start, end, sep = "-")) %>% 
  select(c(20, 18,19, 16, 17, 14, 15, 12, 13)) %>% 
  column_to_rownames("PeakID")

cpm_ChIP_JUNB <- cpm(cts)
cpm_ChIP_JUNB <- as.data.frame(cpm_ChIP_JUNB)

log2cpm_ChIP_JUNB <- log2(cpm_ChIP_JUNB + 1)
log2cpm_ChIP_JUNB <- as.data.frame(log2cpm_ChIP_JUNB)

cpm_ChIP_JUNB <- rownames_to_column(.data = cpm_ChIP_JUNB, var = "peaks")
# write.csv(x = cpm_ChIP_JUNB, file = "data/cpm.csv", row.names = T)

log2cpm_ChIP_JUNB <- rownames_to_column(.data = log2cpm_ChIP_JUNB, var = "peaks")
# write.csv(x = log2cpm_ChIP_JUNB, file = "data/log2cpm.csv", row.names = T)

#groups

# Define the conditions and their replicates (3 replicates each)
replicates <- 2

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
                 "WT_CCM_vs_WT_VEH" = "WT_CCM - WT_VEH")

contrast_numbers <- list(KD_VEH_vs_WT_VEH = 1,
                         KD_VEH_vs_WT_CCM = 2, 
                         KD_CCM_vs_KD_VEH = 3,
                         KD_CCM_vs_WT_VEH = 4,
                         KD_CCM_vs_WT_CCM = 5,
                         WT_CCM_vs_WT_VEH = 6)

cts_matches <- c("KD_V1|KD_V2|WT_V1|WT_V2",  # KD_VEH_vs_WT_VEH
                 "KD_V1|KD_V2|WT_C1|WT_C2",  # KD_VEH_vs_WT_CCM
                 "KD_C1|KD_C2|KD_V1|KD_V2",  # KD_CCM_vs_KD_VEH
                 "KD_C1|KD_C2|WT_V1|WT_V2",  # KD_CCM_vs_WT_VEH
                 "KD_C1|KD_C2|WT_C1|WT_C2",  # KD_CCM_vs_WT_CCM
                 "WT_C1|WT_C2|WT_V1|WT_V2")   # WT_CCT_vs_WT_VEH 

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
  
  cpm_ChIP_JUNB <- apply(cts_tmp, 2, function(x) {(x/sum(x)*1E6)})
  cpm_ChIP_JUNB <- as.data.frame(cpm_ChIP_JUNB)
  
  log2cpm_ChIP_JUNB <- log2(cpm_ChIP_JUNB + 1)
  log2cpm_ChIP_JUNB <- as.data.frame(log2cpm_ChIP_JUNB)
  
  cpm_ChIP_JUNB <- rownames_to_column(.data = cpm_ChIP_JUNB, var = "peaks") %>% 
    filter(peaks %in% peaks_comp)
  
  log2cpm_ChIP_JUNB <- rownames_to_column(.data = log2cpm_ChIP_JUNB, var = "peaks") %>% 
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
  
  
  peak.contrasts.l[[contrast_name]] <- fit_contrast
  peak.contrasts.l[[contrast_name]]$comparison <- paste0(contrast_name)
  
  print(paste0("working on extracting results ", loop_var))
  # Extract top tags
  results_table <- topTags(fit_contrast, adjust.method = "fdr", n = nrow(d1))$table 
  results_table <- results_table %>% 
    rename_all(~ paste0(., "_", contrast_name)) %>% rownames_to_column("peaks") %>% 
    left_join(log2cpm_ChIP_JUNB, by = "peaks")
  
  
  peak.dge.celltype.l[[contrast_name]] <- results_table
  
  # Save results to a CSV file
  file_name <- paste0(deg_folder, "/", contrast_name, ".csv")
  write.csv(results_table, file = file_name, row.names = FALSE)
  # 
}


# peak.dge.celltype.df <- purrr::reduce(peak.dge.celltype.l, full_join, by = "peaks")
# write.csv(peak.dge.celltype.df, file = "../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/deg_tables/all_contrasts_ChIP_JUNB.csv", row.names = FALSE)



dge_values <- names(peak.dge.celltype.l)

for (loop_var in dge_values) {
  data <- peak.dge.celltype.l[[loop_var]]
  data_diff <- data %>% 
    mutate(diffbound = ifelse(data[ ,6] <= 0.05 & data[ ,2] >= 1, "Up", 
                                  ifelse(data[ ,6] <= 0.05 & data[ ,2] <= -1, "Down", "No"))) %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var))))
  
  
  data_diff_filt <- data_diff %>% 
    dplyr::filter(!diffbound == "No")
  
  write.csv(x = data_diff_filt, file = sprintf("../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/deg_difexp_tables/%s_diffbound.csv", loop_var), row.names = F)
  write.csv(x = data_diff, file = sprintf("../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/deg_difexp_tables/%s_diffbound_full.csv", loop_var), row.names = F)
  
}

### ChIP_H3K27ac ----


# Create a folder to save the CSV files

deg_folder <- "../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/deg_tables"
if (!dir.exists(deg_folder)){
  dir.create(deg_folder, recursive = TRUE)
}

deg_difexp_folder <- "../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/deg_difexp_tables"
if (!dir.exists(deg_difexp_folder)){
  dir.create(deg_difexp_folder, recursive = TRUE)
}


# graph volcano folder
if (!dir.exists("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/graphs_ChIP_H3K27ac/volcano")){
  dir.create("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/graphs_ChIP_H3K27ac/volcano", recursive = T)
}


###union peaksets
union_peakset_H3K27ac <- read.table("../results/Human_cells/ChipSeq/data/H3K27ac_mergedpeaks/H3K27ac.filteredNfixed.union.peakSet", sep = "\t")

colnames(union_peakset_H3K27ac) <- c("seqnames",	"start", "end",	"width",	"strand",	"score",	"name",	"label",	"GC",	"N",	"spm")

union_peakset_H3K27ac <- union_peakset_H3K27ac %>% 
  filter(spm > 2) %>% 
  mutate("PeakID" = paste(seqnames, start, end, sep = "-"))

peaks_tofilt <- union_peakset_H3K27ac$PeakID


files <- list.files("../results/Human_cells/ChipSeq/data/H3K27ac_mergedpeaks/peakset_by_comparisons/")
samples <- gsub(".filteredNfixed.union.peakset", "", files[c(5, 4, 1, 3, 2, 6)])


comparison_names <- c("KD_VEH_vs_WT_VEH", "KD_VEH_vs_WT_CCM",
                      "KD_CCM_vs_KD_VEH", "KD_CCM_vs_WT_VEH", 
                      "KD_CCM_vs_WT_CCM", "WT_CCM_vs_WT_VEH")

names(comparison_names) <- samples

peaks_by_comp <- list()
for (sample in samples) {
  
  comparison_name <- comparison_names[[sample]]
  
  tmp <-  read.table(sprintf("../results/Human_cells/ChipSeq/data/H3K27ac_mergedpeaks/peakset_by_comparisons/%s.filteredNfixed.union.peakset", sample), sep = "\t", header = TRUE) %>% 
    mutate("PeakID" = paste(chr, start, end, sep = "-")) %>%  
    filter(PeakID %in% peaks_tofilt) %>% 
    pull("PeakID")
  
  peaks_by_comp[[comparison_name]] <- tmp
  
  rm(tmp, sample, comparison_name)
}


cts <- read.table(file = "../results/Human_cells/ChipSeq/data/multicov_H3K27ac.tsv", sep = "\t")

colnames(cts) <- c("seqnames",	"start", "end",	"width",	"strand",	"score",	"name",	"label",	"GC",	"N",	"spm", 
                   "KD_C1", "KD_C2", "KD_V1", "KD_V2", "WT_C1", "WT_C2", "WT_V1", "WT_V2")

cts <-  cts %>% 
  mutate("PeakID" = paste(seqnames, start, end, sep = "-")) %>% 
  select(c(20, 18,19, 16, 17, 14, 15, 12, 13)) %>% 
  column_to_rownames("PeakID")

cpm_ChIP_H3K27ac <- cpm(cts)
cpm_ChIP_H3K27ac <- as.data.frame(cpm_ChIP_H3K27ac)

log2cpm_ChIP_H3K27ac <- log2(cpm_ChIP_H3K27ac + 1)
log2cpm_ChIP_H3K27ac <- as.data.frame(log2cpm_ChIP_H3K27ac)

cpm_ChIP_H3K27ac <- rownames_to_column(.data = cpm_ChIP_H3K27ac, var = "peaks")
# write.csv(x = cpm_ChIP_H3K27ac, file = "data/cpm.csv", row.names = T)

log2cpm_ChIP_H3K27ac <- rownames_to_column(.data = log2cpm_ChIP_H3K27ac, var = "peaks")
# write.csv(x = log2cpm_ChIP_H3K27ac, file = "data/log2cpm.csv", row.names = T)

#groups

# Define the conditions and their replicates (3 replicates each)
replicates <- 2

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
                 "WT_CCM_vs_WT_VEH" = "WT_CCM - WT_VEH")

contrast_numbers <- list(KD_VEH_vs_WT_VEH = 1,
                         KD_VEH_vs_WT_CCM = 2, 
                         KD_CCM_vs_KD_VEH = 3,
                         KD_CCM_vs_WT_VEH = 4,
                         KD_CCM_vs_WT_CCM = 5,
                         WT_CCM_vs_WT_VEH = 6)

cts_matches <- c("KD_V1|KD_V2|WT_V1|WT_V2",  # KD_VEH_vs_WT_VEH
                 "KD_V1|KD_V2|WT_C1|WT_C2",  # KD_VEH_vs_WT_CCM
                 "KD_C1|KD_C2|KD_V1|KD_V2",  # KD_CCM_vs_KD_VEH
                 "KD_C1|KD_C2|WT_V1|WT_V2",  # KD_CCM_vs_WT_VEH
                 "KD_C1|KD_C2|WT_C1|WT_C2",  # KD_CCM_vs_WT_CCM
                 "WT_C1|WT_C2|WT_V1|WT_V2")   # WT_CCT_vs_WT_VEH 

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
  
  cpm_ChIP_H3K27ac <- apply(cts_tmp, 2, function(x) {(x/sum(x)*1E6)})
  cpm_ChIP_H3K27ac <- as.data.frame(cpm_ChIP_H3K27ac)
  
  log2cpm_ChIP_H3K27ac <- log2(cpm_ChIP_H3K27ac + 1)
  log2cpm_ChIP_H3K27ac <- as.data.frame(log2cpm_ChIP_H3K27ac)
  
  cpm_ChIP_H3K27ac <- rownames_to_column(.data = cpm_ChIP_H3K27ac, var = "peaks") %>% 
    filter(peaks %in% peaks_comp)

  log2cpm_ChIP_H3K27ac <- rownames_to_column(.data = log2cpm_ChIP_H3K27ac, var = "peaks") %>% 
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
  
  
  peak.contrasts.l[[contrast_name]] <- fit_contrast
  peak.contrasts.l[[contrast_name]]$comparison <- paste0(contrast_name)
  
  print(paste0("working on extracting results ", loop_var))
  # Extract top tags
  results_table <- topTags(fit_contrast, adjust.method = "fdr", n = nrow(d1))$table 
  results_table <- results_table %>% 
    rename_all(~ paste0(., "_", contrast_name)) %>% rownames_to_column("peaks") %>% 
    left_join(log2cpm_ChIP_H3K27ac, by = "peaks")
  
  
  peak.dge.celltype.l[[contrast_name]] <- results_table
  
  # Save results to a CSV file
  file_name <- paste0(deg_folder, "/", contrast_name, ".csv")
  write.csv(results_table, file = file_name, row.names = FALSE)
  # 
}


# peak.dge.celltype.df <- purrr::reduce(peak.dge.celltype.l, full_join, by = "peaks")
# write.csv(peak.dge.celltype.df, file = "../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/deg_tables/all_contrasts_ChIP_H3K27ac.csv", row.names = FALSE)



dge_values <- names(peak.dge.celltype.l)

for (loop_var in dge_values) {
  data <- peak.dge.celltype.l[[loop_var]]
  data_diff <- data %>% 
    mutate(diffbound = ifelse(data[ ,6] <= 0.05 & data[ ,2] >= 1, "Up", 
                                  ifelse(data[ ,6] <= 0.05 & data[ ,2] <= -1, "Down", "No"))) %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var))))
  
  
  data_diff_filt <- data_diff %>% 
    dplyr::filter(!diffbound == "No")
  
  write.csv(x = data_diff_filt, file = sprintf("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/deg_difexp_tables/%s_diffbound.csv", loop_var), row.names = F)
  write.csv(x = data_diff, file = sprintf("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/deg_difexp_tables/%s_diffbound_full.csv", loop_var), row.names = F)
  
}

### ChIP_H3K4me1 ----


# Create a folder to save the CSV files

deg_folder <- "../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/deg_tables"
if (!dir.exists(deg_folder)){
  dir.create(deg_folder, recursive = TRUE)
}

deg_difexp_folder <- "../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/deg_difexp_tables"
if (!dir.exists(deg_difexp_folder)){
  dir.create(deg_difexp_folder, recursive = TRUE)
}



# graph volcano folder
if (!dir.exists("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/graphs_ChIP_H3K4me1/volcano")){
  dir.create("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/graphs_ChIP_H3K4me1/volcano", recursive = T)
}


###union peaksets
union_peakset_H3K4me1 <- read.table("../results/Human_cells/ChipSeq/data/H3K4me1_mergedpeaks/H3K4me1.filteredNfixed.union.peakSet", sep = "\t")

colnames(union_peakset_H3K4me1) <- c("seqnames",	"start", "end",	"width",	"strand",	"score",	"name",	"label",	"GC",	"N",	"spm")

union_peakset_H3K4me1 <- union_peakset_H3K4me1 %>% 
  filter(spm > 2) %>% 
  mutate("PeakID" = paste(seqnames, start, end, sep = "-"))

peaks_tofilt <- union_peakset_H3K4me1$PeakID


files <- list.files("../results/Human_cells/ChipSeq/data/H3K4me1_mergedpeaks/peakset_by_comparisons/")
samples <- gsub(".filteredNfixed.union.peakset", "", files[c(5, 4, 1, 3, 2, 6)])


comparison_names <- c("KD_VEH_vs_WT_VEH", "KD_VEH_vs_WT_CCM",
                      "KD_CCM_vs_KD_VEH", "KD_CCM_vs_WT_VEH", 
                      "KD_CCM_vs_WT_CCM", "WT_CCM_vs_WT_VEH")

names(comparison_names) <- samples

peaks_by_comp <- list()
for (sample in samples) {
  
  comparison_name <- comparison_names[[sample]]
  
  tmp <-  read.table(sprintf("../results/Human_cells/ChipSeq/data/H3K4me1_mergedpeaks/peakset_by_comparisons/%s.filteredNfixed.union.peakset", sample), sep = "\t", header = TRUE) %>% 
    mutate("PeakID" = paste(chr, start, end, sep = "-")) %>%  
    filter(PeakID %in% peaks_tofilt) %>% 
    pull("PeakID")
  
  peaks_by_comp[[comparison_name]] <- tmp
  
  rm(tmp, sample, comparison_name)
}


cts <- read.table(file = "../results/Human_cells/ChipSeq/data/multicov_H3K4me1.tsv", sep = "\t")

colnames(cts) <- c("seqnames",	"start", "end",	"width",	"strand",	"score",	"name",	"label",	"GC",	"N",	"spm", 
                   "KD_C1", "KD_C2", "KD_V1", "KD_V2", "WT_C1", "WT_C2", "WT_V1", "WT_V2")

cts <-  cts %>% 
  mutate("PeakID" = paste(seqnames, start, end, sep = "-")) %>% 
  select(c(20, 18,19, 16, 17, 14, 15, 12, 13)) %>% 
  column_to_rownames("PeakID")

cpm_ChIP_H3K4me1 <- cpm(cts)
cpm_ChIP_H3K4me1 <- as.data.frame(cpm_ChIP_H3K4me1)

log2cpm_ChIP_H3K4me1 <- log2(cpm_ChIP_H3K4me1 + 1)
log2cpm_ChIP_H3K4me1 <- as.data.frame(log2cpm_ChIP_H3K4me1)

cpm_ChIP_H3K4me1 <- rownames_to_column(.data = cpm_ChIP_H3K4me1, var = "peaks")
# write.csv(x = cpm_ChIP_H3K4me1, file = "data/cpm.csv", row.names = T)

log2cpm_ChIP_H3K4me1 <- rownames_to_column(.data = log2cpm_ChIP_H3K4me1, var = "peaks")
# write.csv(x = log2cpm_ChIP_H3K4me1, file = "data/log2cpm.csv", row.names = T)


#groups

# Define the conditions and their replicates (3 replicates each)
replicates <- 2

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
                 "WT_CCM_vs_WT_VEH" = "WT_CCM - WT_VEH")

contrast_numbers <- list(KD_VEH_vs_WT_VEH = 1,
                         KD_VEH_vs_WT_CCM = 2, 
                         KD_CCM_vs_KD_VEH = 3,
                         KD_CCM_vs_WT_VEH = 4,
                         KD_CCM_vs_WT_CCM = 5,
                         WT_CCM_vs_WT_VEH = 6)

cts_matches <- c("KD_V1|KD_V2|WT_V1|WT_V2",  # KD_VEH_vs_WT_VEH
                 "KD_V1|KD_V2|WT_C1|WT_C2",  # KD_VEH_vs_WT_CCM
                 "KD_C1|KD_C2|KD_V1|KD_V2",  # KD_CCM_vs_KD_VEH
                 "KD_C1|KD_C2|WT_V1|WT_V2",  # KD_CCM_vs_WT_VEH
                 "KD_C1|KD_C2|WT_C1|WT_C2",  # KD_CCM_vs_WT_CCM
                 "WT_C1|WT_C2|WT_V1|WT_V2")   # WT_CCT_vs_WT_VEH 

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
  
  cpm_ChIP_H3K4me1 <- apply(cts_tmp, 2, function(x) {(x/sum(x)*1E6)})
  cpm_ChIP_H3K4me1 <- as.data.frame(cpm_ChIP_H3K4me1)
  
  log2cpm_ChIP_H3K4me1 <- log2(cpm_ChIP_H3K4me1 + 1)
  log2cpm_ChIP_H3K4me1 <- as.data.frame(log2cpm_ChIP_H3K4me1)
  
  cpm_ChIP_H3K4me1 <- rownames_to_column(.data = cpm_ChIP_H3K4me1, var = "peaks") %>% 
    filter(peaks %in% peaks_comp)
  
 
  log2cpm_ChIP_H3K4me1 <- rownames_to_column(.data = log2cpm_ChIP_H3K4me1, var = "peaks") %>% 
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
  
  
  peak.contrasts.l[[contrast_name]] <- fit_contrast
  peak.contrasts.l[[contrast_name]]$comparison <- paste0(contrast_name)
  
  print(paste0("working on extracting results ", loop_var))
  # Extract top tags
  results_table <- topTags(fit_contrast, adjust.method = "fdr", n = nrow(d1))$table 
  results_table <- results_table %>% 
    rename_all(~ paste0(., "_", contrast_name)) %>% rownames_to_column("peaks") %>% 
    left_join(log2cpm_ChIP_H3K4me1, by = "peaks")
  
  
  peak.dge.celltype.l[[contrast_name]] <- results_table
  
  # Save results to a CSV file
  file_name <- paste0(deg_folder, "/", contrast_name, ".csv")
  write.csv(results_table, file = file_name, row.names = FALSE)
  # 
}


# peak.dge.celltype.df <- purrr::reduce(peak.dge.celltype.l, full_join, by = "peaks")
# write.csv(peak.dge.celltype.df, file = "../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/deg_tables/all_contrasts_ChIP_H3K4me1.csv", row.names = FALSE)



dge_values <- names(peak.dge.celltype.l)

for (loop_var in dge_values) {
  data <- peak.dge.celltype.l[[loop_var]]
  data_diff <- data %>% 
    mutate(diffbound = ifelse(data[ ,6] <= 0.05 & data[ ,2] >= 1, "Up", 
                                  ifelse(data[ ,6] <= 0.05 & data[ ,2] <= -1, "Down", "No"))) %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var))))
  
  
  data_diff_filt <- data_diff %>% 
    dplyr::filter(!diffbound == "No")
  
  write.csv(x = data_diff_filt, file = sprintf("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/deg_difexp_tables/%s_diffbound.csv", loop_var), row.names = F)
  write.csv(x = data_diff, file = sprintf("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/deg_difexp_tables/%s_diffbound_full.csv", loop_var), row.names = F)
  
}


# Optionally save session info
sink("../logs/sessioninfo_human_13_DB_CHIP_edgeR.txt")
sessionInfo()
sink()

