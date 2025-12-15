pacman::p_load( "limma", "tibble", "viridis", "tidyr", "readr", "EnhancedVolcano", "enrichR",
                "RColorBrewer", "pheatmap", "circlize", "stringr", "ComplexHeatmap",
                "dplyr", "ggplot2", "textshape", "edgeR", "paletteer", "forcats")

### ATAC HPAEC ----

# Create a folder to save the CSV files

deg_folder <- "../results/Human_cells/HPAEC_ATAC/edgeR/deg_tables"
deg_difexp_folder <- "../results/Human_cells/HPAEC_ATAC/edgeR/deg_difexp_tables"
graphs <- "../results/Human_cells/HPAEC_ATAC/edgeR/plots/volcano"

# graph volcano folder

dir.create(deg_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_difexp_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(graphs, recursive = TRUE, showWarnings = FALSE)

comparison_names <- c("HPAEC_CCM_vs_HPAEC_VEH")

multicov_colnames.1 <- colnames(read.table("../results/Human_cells/HPAEC_ATAC/HPAEC_merged_peaks/celltype_condition.filteredNfixed.union.peakSet",
                                           sep = "\t", header = TRUE))[1:3]

multicov_colnames.2 <- c("HPAEC_C1", "HPAEC_C2",
                           "HPAEC_V1", "HPAEC_V2")


multicov_hpaec <- read.table("../results/Human_cells/HPAEC_ATAC/data/multicov_output_HPAEC.tsv", sep = "\t")

colnames(multicov_hpaec) <- c(multicov_colnames.1, multicov_colnames.2)

cts <- multicov_hpaec %>% 
  mutate(peaks = paste(seqnames, start, end, sep = "-"), .keep = "unused") %>% 
  column_to_rownames("peaks") %>% 
  select(starts_with("HPAEC_V"), starts_with("HPAEC_C")) 


write_csv(cts, file = "../results/Human_cells/HPAEC_ATAC/data/cts_HPAEC.csv")

cpm_atac <- cpm(cts)
cpm_atac <- as.data.frame(cpm_atac)

log2cpm_atac <- log2(cpm_atac + 1)
log2cpm_atac <- as.data.frame(log2cpm_atac)

cpm_atac <- rownames_to_column(.data = cpm_atac, var = "peaks")
write.csv(x = cpm_atac, file = "../results/Human_cells/HPAEC_ATAC/data/cpm.csv", row.names = F)


log2cpm_atac <- rownames_to_column(.data = log2cpm_atac, var = "peaks")
write.csv(x = log2cpm_atac, file = "../results/Human_cells/HPAEC_ATAC/data/log2cpm.csv", row.names = F)

#groups

ctsGroups <- c("HPAEC_VEH", "HPAEC_VEH",
               "HPAEC_CCM", "HPAEC_CCM")
  
  
d <- DGEList(counts = cts,
             group = factor(ctsGroups))
  
d

d <- calcNormFactors(d)
  
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
  
d1 <- estimateDisp(d, design.mat)
dim(d1$counts)
  
# Fit the negative binomial GLM
fit <- glmFit(d1, design.mat)


# Perform LRT for the current contrast
lrt <- glmLRT(fit, contrast = makeContrasts(HPAEC_CCM - HPAEC_VEH, levels = design.mat))
  
reslrt <- topTags(lrt, n = Inf, adjust.method = "fdr")$table %>%   # or topTags(qlt, n=Inf)$table
  rownames_to_column("peaks") %>%
  rename_with(~ paste0(., "_HPAEC_ccm_vs_HPAEC_veh"), -peaks)

reslrt_sigs <- reslrt %>%
  mutate(diffaccessible = case_when(FDR_HPAEC_ccm_vs_HPAEC_veh <= 0.05 & logFC_HPAEC_ccm_vs_HPAEC_veh >= 1 ~ "Up",
                                   FDR_HPAEC_ccm_vs_HPAEC_veh <= 0.05 & logFC_HPAEC_ccm_vs_HPAEC_veh <= -1 ~ "Down",
                                   TRUE ~ "No"  )) %>%
  arrange(desc(logFC_HPAEC_ccm_vs_HPAEC_veh))


reslrt_sigs %>% group_by(diffaccessible) %>% summarise(DAR = n())

write_csv(filter(reslrt), file = file.path(deg_folder, "HPAEC_ccm_vs_HPAEC_veh.csv"))
write_csv(filter(reslrt_sigs, diffaccessible != "No"), file = file.path(deg_difexp_folder, "HPAEC_ccm_vs_HPAEC_veh_diffacc.csv"))
write_csv(reslrt_sigs, file = file.path(deg_difexp_folder, "HPAEC_ccm_vs_HPAEC_veh_diffacc_full.csv"))


### ATAC HUVEC ----

# Create a folder to save the CSV files

deg_folder <- "../results/Human_cells/HUVEC_ATAC/edgeR/deg_tables"
deg_difexp_folder <- "../results/Human_cells/HUVEC_ATAC/edgeR/deg_difexp_tables"
graphs <- "../results/Human_cells/HUVEC_ATAC/edgeR/plots/volcano"

# graph volcano folder

dir.create(deg_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_difexp_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(graphs, recursive = TRUE, showWarnings = FALSE)

comparison_names <- c("HUVEC_CCM_vs_HUVEC_VEH")

multicov_colnames.1 <- colnames(read.table("../results/Human_cells/HUVEC_ATAC/HUVEC_merged_peaks/celltype_condition.filteredNfixed.union.peakSet",
                                           sep = "\t", header = TRUE))[1:3]

multicov_colnames.2 <- c("HUVEC_C1", "HUVEC_C2",
                         "HUVEC_V1", "HUVEC_V2")


multicov_HUVEC <- read.table("../results/Human_cells/HUVEC_ATAC/data/multicov_output_HUVEC.tsv", sep = "\t")

colnames(multicov_HUVEC) <- c(multicov_colnames.1, multicov_colnames.2)

cts <- multicov_HUVEC %>% 
  mutate(peaks = paste(seqnames, start, end, sep = "-"), .keep = "unused") %>% 
  column_to_rownames("peaks") %>% 
  select(starts_with("HUVEC_V"), starts_with("HUVEC_C")) 


write_csv(cts, file = "../results/Human_cells/HUVEC_ATAC/data/cts_HUVEC.csv")

cpm_atac <- cpm(cts)
cpm_atac <- as.data.frame(cpm_atac)

log2cpm_atac <- log2(cpm_atac + 1)
log2cpm_atac <- as.data.frame(log2cpm_atac)

cpm_atac <- rownames_to_column(.data = cpm_atac, var = "peaks")
write.csv(x = cpm_atac, file = "../results/Human_cells/HUVEC_ATAC/data/cpm.csv", row.names = F)


log2cpm_atac <- rownames_to_column(.data = log2cpm_atac, var = "peaks")
write.csv(x = log2cpm_atac, file = "../results/Human_cells/HUVEC_ATAC/data/log2cpm.csv", row.names = F)

#groups

ctsGroups <- c("HUVEC_VEH", "HUVEC_VEH",
               "HUVEC_CCM", "HUVEC_CCM")


d <- DGEList(counts = cts,
             group = factor(ctsGroups))

d

d <- calcNormFactors(d)

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)

d1 <- estimateDisp(d, design.mat)
dim(d1$counts)

# Fit the negative binomial GLM
fit <- glmFit(d1, design.mat)


# Perform LRT for the current contrast
lrt <- glmLRT(fit, contrast = makeContrasts(HUVEC_CCM - HUVEC_VEH, levels = design.mat))

reslrt <- topTags(lrt, n = Inf, adjust.method = "fdr")$table %>%   # or topTags(qlt, n=Inf)$table
  rownames_to_column("peaks") %>%
  rename_with(~ paste0(., "_HUVEC_ccm_vs_HUVEC_veh"), -peaks)

reslrt_sigs <- reslrt %>%
  mutate(diffaccessible = case_when(FDR_HUVEC_ccm_vs_HUVEC_veh <= 0.05 & logFC_HUVEC_ccm_vs_HUVEC_veh >= 1 ~ "Up",
                                    FDR_HUVEC_ccm_vs_HUVEC_veh <= 0.05 & logFC_HUVEC_ccm_vs_HUVEC_veh <= -1 ~ "Down",
                                    TRUE ~ "No"  )) %>%
  arrange(desc(logFC_HUVEC_ccm_vs_HUVEC_veh))


reslrt_sigs %>% group_by(diffaccessible) %>% summarise(DAR = n())

write_csv(filter(reslrt), file = file.path(deg_folder, "HUVEC_ccm_vs_HUVEC_veh.csv"))
write_csv(filter(reslrt_sigs, diffaccessible != "No"), file = file.path(deg_difexp_folder, "HUVEC_ccm_vs_HUVEC_veh_diffacc.csv"))
write_csv(reslrt_sigs, file = file.path(deg_difexp_folder, "HUVEC_ccm_vs_HUVEC_veh_diffacc_full.csv"))


### RNA HPAEC ----

# Create a folder to save the CSV files

deg_folder_rna <- "../results/Human_cells/HPAEC_RNA/edgeR/deg_tables"
deg_difexp_folder_rna <- "../results/Human_cells/HPAEC_RNA/edgeR/deg_difexp_tables"
graphs <- "../results/Human_cells/HPAEC_RNA/edgeR/plots/volcano"
data_folder_rna <- "../results/Human_cells/HPAEC_RNA/data"
# graph volcano folder

dir.create(deg_folder_rna, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_difexp_folder_rna, recursive = TRUE, showWarnings = FALSE)
dir.create(graphs, recursive = TRUE, showWarnings = FALSE)
dir.create(data_folder_rna, recursive = TRUE, showWarnings = FALSE)

comparison_names <- c("HPAEC_CCM_vs_HPAEC_VEH")

counts <- read.table("../results/Human_cells/data/HPAEC_HUVEC_counts_length_scaled.tsv", header = TRUE)

colnames(counts)

cts <- counts %>% 
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames("gene_name") %>% 
  select(starts_with("HA_V"), starts_with("HA_C")) 


write.csv(cts, file = "../results/Human_cells/HPAEC_RNA/data/cts_HPAEC.csv")

cpm_rna <- cpm(cts)
cpm_rna <- as.data.frame(cpm_rna)

log2cpm_rna <- log2(cpm_rna + 1)
log2cpm_rna <- as.data.frame(log2cpm_rna)

cpm_rna <- rownames_to_column(.data = cpm_rna, var = "gene")
write.csv(x = cpm_rna, file = "../results/Human_cells/HPAEC_RNA/data/cpm.csv", row.names = F)

log2cpm_rna <- rownames_to_column(.data = log2cpm_rna, var = "gene")
write.csv(x = log2cpm_rna, file = "../results/Human_cells/HPAEC_RNA/data/log2cpm.csv", row.names = F)

#groups

ctsGroups <- c("HPAEC_VEH", "HPAEC_VEH", "HPAEC_VEH",
               "HPAEC_CCM", "HPAEC_CCM", "HPAEC_CCM")


d <- DGEList(counts = cts,
             group = factor(ctsGroups))


dim(d)

# avg_logCPM <- rowMeans(log2cpm_rna[ ,-1])
# keep <- avg_logCPM > 1.5 
keep <- filterByExpr(d)

cpm_rna_f <- cpm_rna[keep,]
log2cpm_rna_f <- log2cpm_rna[keep,]

write.csv(x = cpm_rna_f, file = "../results/Human_cells/HPAEC_RNA/data/cpm_filt.csv", row.names = F)
write.csv(x = log2cpm_rna_f, file = "../results/Human_cells/HPAEC_RNA/data/log2cpm_filt.csv", row.names = F)

d <- d[keep,]
dim(d)

d <- calcNormFactors(d)

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)

d1 <- estimateDisp(d, design.mat)
dim(d1$counts)

# Fit the negative binomial GLM
fit <- glmFit(d1, design.mat)


# Perform LRT for the current contrast
lrt <- glmLRT(fit, contrast = makeContrasts(HPAEC_CCM - HPAEC_VEH, levels = design.mat))

reslrt <- topTags(lrt, n = Inf, adjust.method = "fdr")$table %>%   # or topTags(qlt, n=Inf)$table
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_HPAEC_ccm_vs_HPAEC_veh"), -gene)

reslrt_sigs <- reslrt %>%
  mutate(diffexpressed = case_when(FDR_HPAEC_ccm_vs_HPAEC_veh <= 0.05 & logFC_HPAEC_ccm_vs_HPAEC_veh >= 0.5 ~ "Up",
                                   FDR_HPAEC_ccm_vs_HPAEC_veh <= 0.05 & logFC_HPAEC_ccm_vs_HPAEC_veh <= -0.5 ~ "Down",
                                   TRUE ~ "No"  )) %>%
  arrange(desc(logFC_HPAEC_ccm_vs_HPAEC_veh))


reslrt_sigs %>% group_by(diffexpressed) %>% summarise(DEG = n())

write_csv(filter(reslrt), file = file.path(deg_folder_rna, "HPAEC_ccm_vs_HPAEC_veh.csv"))
write_csv(filter(reslrt_sigs, diffexpressed != "No") %>% left_join(cts %>% rownames_to_column("gene"), by = "gene"), file = file.path(deg_difexp_folder_rna, "HPAEC_ccm_vs_HPAEC_veh_diffexp.csv"))
write_csv(reslrt_sigs %>% left_join(cts %>% rownames_to_column("gene"), by = "gene"), file = file.path(deg_difexp_folder_rna, "HPAEC_ccm_vs_HPAEC_veh_diffexp_full.csv"))


### RNA HUVEC ----

# Create a folder to save the CSV files

deg_folder_rna <- "../results/Human_cells/HUVEC_RNA/edgeR/deg_tables"
deg_difexp_folder_rna <- "../results/Human_cells/HUVEC_RNA/edgeR/deg_difexp_tables"
graphs <- "../results/Human_cells/HUVEC_RNA/edgeR/plots/volcano"
data_folder_rna <- "../results/Human_cells/HUVEC_RNA/data"
# graph volcano folder

dir.create(deg_folder_rna, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_difexp_folder_rna, recursive = TRUE, showWarnings = FALSE)
dir.create(graphs, recursive = TRUE, showWarnings = FALSE)
dir.create(data_folder_rna, recursive = TRUE, showWarnings = FALSE)

cts <-  counts %>% 
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames("gene_name") %>% 
  select(starts_with("HV_V"), starts_with("HV_C")) 



write.csv(cts, file = "../results/Human_cells/HUVEC_RNA/data/cts_HUVEC.csv")

cpm_rna <- cpm(cts)
cpm_rna <- as.data.frame(cpm_rna)

log2cpm_rna <- log2(cpm_rna + 1)
log2cpm_rna <- as.data.frame(log2cpm_rna)

cpm_rna <- rownames_to_column(.data = cpm_rna, var = "gene")
write.csv(x = cpm_rna, file = "../results/Human_cells/HUVEC_RNA/data/cpm.csv", row.names = F)


log2cpm_rna <- rownames_to_column(.data = log2cpm_rna, var = "gene")
write.csv(x = log2cpm_rna, file = "../results/Human_cells/HUVEC_RNA/data/log2cpm.csv", row.names = F)

#groups

ctsGroups <- c("HUVEC_VEH", "HUVEC_VEH", "HUVEC_VEH", 
               "HUVEC_CCM", "HUVEC_CCM", "HUVEC_CCM")


d <- DGEList(counts = cts,
             group = factor(ctsGroups))

dim(d)
# avg_logCPM <- rowMeans(log2cpm_rna[ ,-1])
# keep <- avg_logCPM > 1.5  # Adjust this based on your histogram
keep <- filterByExpr(d)
cpm_rna_f <- cpm_rna[keep,]
log2cpm_rna_f <- log2cpm_rna[keep,]

write.csv(x = cpm_rna_f, file = "../results/Human_cells/HUVEC_RNA/data/cpm_filt.csv", row.names = F)
write.csv(x = log2cpm_rna_f, file = "../results/Human_cells/HUVEC_RNA/data/log2cpm_filt.csv", row.names = F)

d <- d[keep,]
dim(d)

d <- calcNormFactors(d)

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)

d1 <- estimateDisp(d, design.mat)
dim(d1$counts)

# Fit the negative binomial GLM
fit <- glmFit(d1, design.mat)


# Perform LRT for the current contrast
lrt <- glmLRT(fit, contrast = makeContrasts(HUVEC_CCM - HUVEC_VEH, levels = design.mat))

reslrt <- topTags(lrt, n = Inf, adjust.method = "fdr")$table %>%   # or topTags(qlt, n=Inf)$table
  rownames_to_column("gene") %>%
  rename_with(~ paste0(., "_HUVEC_ccm_vs_HUVEC_veh"), -gene)

reslrt_sigs <- reslrt %>%
  mutate(diffexpressed = case_when(FDR_HUVEC_ccm_vs_HUVEC_veh <= 0.05 & logFC_HUVEC_ccm_vs_HUVEC_veh >= 0.5 ~ "Up",
                                   FDR_HUVEC_ccm_vs_HUVEC_veh <= 0.05 & logFC_HUVEC_ccm_vs_HUVEC_veh <= -0.5 ~ "Down",
                                   TRUE ~ "No"  )) %>%
  arrange(desc(logFC_HUVEC_ccm_vs_HUVEC_veh))


reslrt_sigs %>% group_by(diffexpressed) %>% summarise(DAR = n())

write_csv(filter(reslrt), file = file.path(deg_folder_rna, "HUVEC_ccm_vs_HUVEC_veh.csv"))
write_csv(filter(reslrt_sigs, diffexpressed != "No") %>% left_join(cts %>% rownames_to_column("gene"), by = "gene"), file = file.path(deg_difexp_folder_rna, "HUVEC_ccm_vs_HUVEC_veh_diffexp.csv"))
write_csv(reslrt_sigs%>% left_join(cts %>% rownames_to_column("gene"), by = "gene"), file = file.path(deg_difexp_folder_rna, "HUVEC_ccm_vs_HUVEC_veh_diffexp_full.csv"))


## volcano ----
hpaec_deg <- read.csv("../results/Human_cells/HPAEC_RNA/edgeR/deg_difexp_tables/HPAEC_ccm_vs_HPAEC_veh_diffexp_full.csv")
huvec_deg <- read.csv("../results/Human_cells/HUVEC_RNA/edgeR/deg_difexp_tables/HUVEC_ccm_vs_HUVEC_veh_diffexp_full.csv")

summary_data_hpaec <- hpaec_deg %>% 
  group_by(diffexpressed) %>% 
  summarise(n = n())
summary_data_huvec <- huvec_deg %>% 
  group_by(diffexpressed) %>% 
  summarise(n = n())

hpaec_volcano_p <- EnhancedVolcano(hpaec_deg, 
                                   lab = NA,
                                   legendLabSize = 12,
                                   axisLabSize = 13, 
                                   x = "logFC_HPAEC_ccm_vs_HPAEC_veh",
                                   y = "FDR_HPAEC_ccm_vs_HPAEC_veh", 
                                   title = "HPAEC\nCCM-like env vs Vehicle",
                                   subtitleLabSize = 17,
                                   subtitle = paste0(summary_data_hpaec[1,]$n, " downregulated and ", summary_data_hpaec[3,]$n, " upregulated"),
                                   pCutoff = 0.05,
                                   FCcutoff = 0.5, 
                                   raster = T) 

huvec_volcano_p <- EnhancedVolcano(huvec_deg, 
                                   lab = NA,
                                   legendLabSize = 12,
                                   axisLabSize = 13, 
                                   x = "logFC_HUVEC_ccm_vs_HUVEC_veh",
                                   y = "FDR_HUVEC_ccm_vs_HUVEC_veh", 
                                   title = "HUVEC\nCCM-like env vs Vehicle",
                                   subtitleLabSize = 17,
                                   subtitle = paste0(summary_data_huvec[1,]$n, " downregulated and ", summary_data_huvec[3,]$n, " upregulated"),
                                   pCutoff = 0.05,
                                   FCcutoff = 0.5, 
                                   raster = T) 

ggsave(filename = "../results/Human_cells/plots/hpaec_volcano.png",
       plot = hpaec_volcano_p, width = 6, height = 6, units = "in", dpi = 300, device = "png")
ggsave(filename = "../results/Human_cells/plots/huvec_volcano.png",
       plot = huvec_volcano_p, width = 6, height = 6, units = "in", dpi = 300, device = "png")

## venn----
hpaec_deg <- read.csv("../results/Human_cells/HPAEC_RNA/edgeR/deg_difexp_tables/HPAEC_ccm_vs_HPAEC_veh_diffexp_full.csv")
huvec_deg <- read.csv("../results/Human_cells/HUVEC_RNA/edgeR/deg_difexp_tables/HUVEC_ccm_vs_HUVEC_veh_diffexp_full.csv")
hcmec_deg <- read.csv("../results/Human_cells/D3_RNA/edgeR/deg_difexp_tables/WT_CCM_vs_WT_VEH_diffexp.csv")

#up and down genes from each DF
genes_up_l <- lapply(list(hpaec = hpaec_deg, huvec = huvec_deg, hcmec = hcmec_deg), function(x){
  x %>%
    filter(diffexpressed == "Up") %>%
    pull(gene)
})

genes_down_l <- lapply(list(hpaec = hpaec_deg, huvec = huvec_deg, hcmec = hcmec_deg), function(x){
  x %>%
    filter(diffexpressed == "Down") %>%
    pull(gene)
})


# 
# from_list <- function(list_data) {
#   members = unique(unlist(list_data))
#   data.frame(
#     lapply(list_data, function(set) members %in% set),
#     row.names=members,
#     check.names=FALSE
#   )
# }
# 
# new_lists_up_matrix <-from_list(genes_up_l)

names(genes_up_l) <- c("HPAEC + CCM-like env",
                       "HUVEC + CCM-like env",
                       "hCMEC + CCM-like env")


# 8332AC
comparison_up_gene <- ggvenn::ggvenn(genes_up_l, fill_color = c("#BAD1CD", "#F2D1C9", "#E086D3"),
                                     set_name_size = 0, text_size = 3.5, show_percentage = FALSE, set_name_color = "white") + 
  ggtitle(label = "Overlapped upregulated genes\nunder CCM-like env") + 
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold", margin = margin(b= 7, t = 2, unit = "pt"))) +
  annotate("text", x = -.85, y = 1.9, label = "HPAEC", size = 3.2) +
  annotate("text", x = 0.85, y = 1.9, label = "HUVEC", size = 3.2) +
  annotate("text", x = 0, y = -1.9, label = "hCMEC", size = 3.2) 

ggsave(filename = "../results/Human_cells/plots/up_gene_list_venn_hcells.png", plot = comparison_up_gene, width = 3.5, height = 3.5, bg = "white")

comparison_down_gene <- ggvenn::ggvenn(genes_down_l, fill_color = c("#BAD1CD", "#F2D1C9", "#E086D3"),
                                       set_name_size = 0, text_size = 3.5, show_percentage = FALSE, set_name_color = "white") + 
  ggtitle(label = "Overlapped downregulated genes\nunder CCM-like env") + 
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold", margin = margin(b= 7, t = 2, unit = "pt"))) +
  annotate("text", x = -.85, y = 1.9, label = "HPAEC", size = 3.2) +
  annotate("text", x = 0.85, y = 1.9, label = "HUVEC", size = 3.2) +
  annotate("text", x = 0, y = -1.9, label = "hCMEC", size = 3.2) 

ggsave(filename = "../results/Human_cells/plots/down_gene_list_venn_hcells.png", plot = comparison_down_gene, width = 3.5, height = 3.5, bg = "white")



#heatmap ----
from_list <- function(list_data) {
  members = unique(unlist(list_data))
  data.frame(
    lapply(list_data, function(set) members %in% set),
    row.names=members,
    check.names=FALSE
  )
}

gene_boolean_df <- from_list(genes_up_l)

##common genes
gene_boolean_df_all <- gene_boolean_df %>%
  filter(if_all(everything(), ~ .x)) %>% 
  rownames_to_column("gene") %>% pull(gene)

# ##Unique per cell type
# # keep only genes that are TRUE in exactly 1 column
# unique_mask <- gene_boolean_df[rowSums(gene_boolean_df) == 1, , drop = FALSE]
# # for each gene, find which column is TRUE
# cell_for_gene <- colnames(unique_mask)[max.col(unique_mask, ties.method = "first")]
# 
# gene_unique_df <- data.frame(gene = rownames(unique_mask),
#                              cell = cell_for_gene,
#                              row.names = NULL,
#                              check.names = FALSE) %>%
#   mutate(cell = dplyr::recode(cell,
#                               "HPAEC + CCM-like env" = "HPAEC",
#                               "HUVEC + CCM-like env" = "HUVEC",
#                               "hCMEC + CCM-like env" = "hCMEC"))
# 
# ### all genes ----
# all_genes <- gene_boolean_df %>% 
#   rownames_to_column("gene") %>% 
#   arrange(desc(`HPAEC + CCM-like env`), desc(`HUVEC + CCM-like env`), desc(`hCMEC + CCM-like env`)) %>% 
#   pull("gene")

hp_cts <- read.csv("../results/Human_cells/HPAEC_RNA/data/log2cpm.csv") %>% filter(gene %in% gene_boolean_df_all)
hv_cts <- read.csv("../results/Human_cells/HUVEC_RNA/data/log2cpm.csv") %>% filter(gene %in% gene_boolean_df_all)
hc_cts <- read.csv("../results/Human_cells/D3_RNA/data/log2cpm.csv") %>% 
  select(gene, contains("WT_V"), contains("WT_C")) %>% filter(gene %in% gene_boolean_df_all)


full_cts <- list(hp_cts, hv_cts, hc_cts) %>% 
  purrr::reduce(full_join, by = "gene") %>% 
  column_to_rownames("gene")


full_mtx <- full_cts 
groups <- c("HA_", "HV_", "WT_")

for (g in groups) {
  cols <- grep(paste0("^", g), colnames(full_cts))
  full_mtx[, cols] <- t(scale(t(full_cts[, cols])))
}

full_mtx[is.na(full_mtx)] <- 0

color_annotmap <- setNames(c("#B2DF8A", "#33A02C"),
                          c("Vehicle\n", "CCM-like\nenvironment"))

top_annot_gene <- HeatmapAnnotation( Condition = rep(c("Vehicle\n", "CCM-like\nenvironment"), each = 3),
                  col = list(Condition = color_annotmap),
                  annotation_height = unit(c(2.5, 2.5), "mm"),  # Adjust heights for bars
                  show_legend = FALSE,
                  show_annotation_name = FALSE)

top_annot_gene2 <- HeatmapAnnotation( Condition = rep(rep(c("Vehicle\n", "CCM-like\nenvironment"), each = 3), 3),
                                     col = list(Condition = color_annotmap),
                                     annotation_height = unit(c(2.5, 2.5), "mm"),  # Adjust heights for bars
                                     show_legend = FALSE,
                                     show_annotation_name = FALSE)

col_fun_rna <- circlize::colorRamp2(seq(-3, 3, length = 100), plasma(100))
hmap_hpaec <- ComplexHeatmap::Heatmap(full_mtx[ , c(1:6)],  
                                          col = col_fun_rna,
                                          show_row_names = FALSE,
                                          show_heatmap_legend = FALSE, 
                                          #column_title = "HPAEC", 
                                          column_title_gp = gpar(fontsize = 7),
                                          show_column_names = FALSE,
                                          cluster_rows = FALSE, 
                                      use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Hanning",
                                          top_annotation = top_annot_gene,
                                          cluster_columns = FALSE)

hmap_huvec <- ComplexHeatmap::Heatmap(full_mtx[ , c(7:12)],  
                                      col = col_fun_rna,
                                      show_row_names = FALSE,
                                      show_heatmap_legend = FALSE, 
                                      #column_title = "HUVEC", 
                                      column_title_gp = gpar(fontsize = 7),
                                      show_column_names = FALSE,
                                      cluster_rows = FALSE, 
                                      use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Hanning",
                                      top_annotation = top_annot_gene,
                                      cluster_columns = FALSE)

hmap_hcmec <- ComplexHeatmap::Heatmap(full_mtx[ , c(13:18)],  
                                      col = col_fun_rna,
                                      show_row_names = F,
                                      show_heatmap_legend = FALSE, 
                                      #column_title = "hCMEC", 
                                      column_title_gp = gpar(fontsize = 7),
                                      show_column_names = FALSE,
                                      cluster_rows = FALSE, 
                                      use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Hanning",
                                      top_annotation = top_annot_gene,
                                      cluster_columns = FALSE)




draw(hmap_hpaec + hmap_huvec + hmap_hcmec)

#### linked genes ----
ccm_matched <- read.csv("../results/Human_cells/timepoint_reg/tables/linked_genes_commonCCM_humanEC.csv")

ccm_matched %>% group_by(classification) %>% summarise(number = n())

hp_cts_enh <- read.csv("../results/Human_cells/HPAEC_RNA/data/log2cpm_filt.csv") %>% filter(gene %in% ccm_matched$linked_gene)
hv_cts_enh <- read.csv("../results/Human_cells/HUVEC_RNA/data/log2cpm_filt.csv") %>% filter(gene %in% ccm_matched$linked_gene)
hc_cts_enh <- read.csv("../results/Human_cells/D3_RNA/data/log2cpm_filt.csv") %>% 
  select(gene, contains("WT_V"), contains("WT_C")) %>% filter(gene %in% ccm_matched$linked_gene)

full_cts_enh <- list(hp_cts_enh, hv_cts_enh, hc_cts_enh) %>% 
  purrr::reduce(full_join, by = "gene") %>% drop_na() %>% 
  column_to_rownames("gene")

full_cts_mtx <- full_cts_enh

groups <- c("HA_", "HV_", "WT_")

for (g in groups) {
  cols <- grep(paste0("^", g), colnames(full_cts_enh))
  full_cts_mtx[, cols] <- t(scale(t(full_cts_enh[, cols])))
}

full_cts_mtx <- full_cts_mtx %>% 
  rownames_to_column("gene")

gene_enh_pairs_hmap <- ccm_matched %>% 
  select(linked_gene, classification) %>% 
  filter(classification == "Common") %>% 
  inner_join(full_cts_mtx, by = join_by("linked_gene" == "gene")) %>%  
  distinct(linked_gene, .keep_all = TRUE) %>% 
  mutate(avg_value = rowMeans(across(starts_with(c("HA_C", "HV_C", "WT_C")))),
         d3avg_value = rowMeans(across(starts_with(c("WT_C"))))) %>% 
  filter(d3avg_value > 0) %>% 
  #mutate(avg_value = rowMeans(across(where(is.numeric)))) %>% 
  arrange(desc(avg_value)) %>% select(-avg_value, -d3avg_value)
# 
# lenghts <- gene_enh_pairs_hmap %>% group_by(classification) %>% summarise(number = n()) %>% deframe()
# 
# partition_hu_enh <- fct_inorder(c(rep("HPAEC",lenghts["HPAEC"]),
#                                   rep("HUVEC", lenghts["HUVEC"]),
#                                   rep("hCMEC/D3", lenghts["hCMEC/D3"]),
#                                   rep("Common", lenghts["Common"])))
# 
# partition_hu_hmap <- Heatmap(partition_hu_enh, col=structure(c("#BAD1CD", "#F2D1C9", "#E086D3",  "#64113F"), 
#                                                                names = c("HPAEC", 
#                                                                          "HUVEC",  
#                                                                          "hCMEC/D3",
#                                                                          "Common")), 
#                            show_column_names = FALSE, name = " ", show_heatmap_legend = FALSE, 
#                            #row_order = order_for_heatmaps_hu,   
#                            row_title_gp = gpar(fontsize = 7),
#                            show_row_names = FALSE, width=unit(1,'mm'))




top_annot_block <- HeatmapAnnotation(Source = anno_block(gp = gpar(fill = c("#B2DF8A", "#33A02C"))),
                                     height = unit(2.5, "mm"),
                                     show_legend = FALSE, 
                                     show_annotation_name = FALSE)



cond_split <- factor(rep(c("Vehicle", "CCM-like env"), each = 3),
                     levels = c("Vehicle", "CCM-like env"))

peak_gene_ids <- gene_enh_pairs_hmap %>%
  dplyr::filter(classification == "Common",
                linked_gene %in% c("SERPINE1", "ICAM1", "CCL2", "VEGFC", "VEGFA", "TGFB2")) %>% 
  mutate(position =  match(linked_gene, gene_enh_pairs_hmap[gene_enh_pairs_hmap$classification == "Common" , "linked_gene"])) %>% 
  arrange(position)



ha <- rowAnnotation(gene_label = anno_mark(at = peak_gene_ids$position,
                                           labels = peak_gene_ids$linked_gene,
                                           labels_gp = gpar(fontsize = 4, col = "black"),
                                           padding = unit(1, "mm"),
                                           link_width = unit(2, "mm")))

col_fun_rna <- circlize::colorRamp2(seq(-2.5, 2.5, length = 100), plasma(100))

hmap_hpaec_enh <- ComplexHeatmap::Heatmap(gene_enh_pairs_hmap[, c(3:8)],  
                                      col = col_fun_rna,
                                      show_row_names = FALSE,
                                      show_heatmap_legend = FALSE, 
                                      column_title = NULL,
                                      show_column_names = FALSE,
                                      cluster_rows = FALSE, border =T,
                                      column_gap = unit(0.5, "mm"),
                                      top_annotation = top_annot_block,
                                      column_split = cond_split,na_col = "white",
                                      cluster_columns = FALSE)

hmap_huvec_enh <- ComplexHeatmap::Heatmap(gene_enh_pairs_hmap[, c(9:14)],  
                                      col = col_fun_rna,
                                      show_row_names = FALSE,
                                      show_heatmap_legend = FALSE, 
                                      column_title = NULL,
                                      show_column_names = FALSE,
                                      cluster_rows = FALSE, na_col = "white",
                                      column_gap = unit(0.5, "mm"),border =T,
                                      top_annotation = top_annot_block,
                                      column_split = cond_split,
                                      cluster_columns = FALSE)

hmap_hcmec_enh <- ComplexHeatmap::Heatmap(gene_enh_pairs_hmap[, c(15:20)],  
                                      col = col_fun_rna,
                                      show_row_names = F, na_col = "white",
                                      show_heatmap_legend = FALSE, 
                                      column_title = NULL, 
                                      show_column_names = FALSE,border =T,
                                      cluster_rows = FALSE, 
                                      column_gap = unit(0.5, "mm"),
                                      top_annotation = top_annot_block,
                                      column_split = cond_split,
                                      right_annotation = ha,
                                      cluster_columns = FALSE)

png("../results/Human_cells/timepoint_reg/heatmaps_coverage/rna_common_linkgenes.png", width = 2.0, height = 1.5, units = "in", res = 300)
draw(hmap_hpaec_enh + hmap_huvec_enh + hmap_hcmec_enh, gap = unit(1, "mm"))
dev.off()


# Set Enrichr site----
setEnrichrSite("Enrichr")

dbs <- listEnrichrDbs()
dbs <- "KEGG_2021_Human"
irrelevant_terms <- read.csv(file = "utils/Irrelevant_terms.csv") %>% pull("Irrelevant_terms")

hp_cts_bg <- read.csv("../results/Human_cells/HPAEC_RNA/data/log2cpm_filt.csv") %>% pull(gene)
hv_cts_bg <- read.csv("../results/Human_cells/HUVEC_RNA/data/log2cpm_filt.csv") %>% pull(gene)
hc_cts_bg <- read.csv("../results/Human_cells/D3_RNA/data/log2cpm_filt.csv") %>% pull(gene)

bg_genes <- unique(c(hp_cts_bg, hv_cts_bg, hc_cts_bg))


enrichr_common <- enrichr(gene_enh_pairs_hmap$linked_gene, dbs, bg_genes, include_overlap = TRUE)

enrich_df <- as.data.frame(enrichr_common$KEGG_2021_Human) %>% 
    filter(P.value <= 0.05,
           !Term %in% irrelevant_terms) %>% 
    separate(col = Overlap, into = c("Count", "total"), sep = "/", convert = T, remove = FALSE) %>% 
    mutate("log10_AdjPval" = -log10(Adjusted.P.value),
           "log10_Pval" = -log10(P.value)) %>%
    select(-contains("Old"))

top_terms <- enrich_df %>% 
  slice_max(order_by = log10_Pval, n = 12) 


enrich_plot <- ggplot(top_terms, aes(x = log10_Pval, y = reorder(Term, log10_Pval), color = log10_Pval, size = Count)) +
  geom_point() +
  scale_color_viridis_c(breaks = seq(floor(min(enrich_df$log10_Pval)), ceiling(max(enrich_df$log10_Pval)), length.out = 5),
                        labels = function(x) round(x, 0),name = "-log10\nP-value") +
  theme_bw() +
  ggtitle("KEGG enrichment") +
  theme(axis.text.y = element_text(size = 9, face = "bold", colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.key.height = unit(12, "pt"),
        legend.title = element_text(size = 10, face = "bold", margin = margin(t = 2, b = 4)),
        legend.text = element_text(size = 9),
        axis.ticks.x = element_line(linewidth = 0),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 9.5, face = "bold", hjust = 0.5, vjust = 0.5, colour = "black")) +
  scale_size_continuous(range = c(2, 4))

ggsave("../results/Human_cells/plots/kegg_linked_genes_common.png", enrich_plot, height = 3.5, width = 4.6, units = "in", dpi = 300)

# Optionally save session info
sink("../logs/sessioninfo_human_19_HPAEC_HUVEC_edger.txt")
sessionInfo()
sink()

