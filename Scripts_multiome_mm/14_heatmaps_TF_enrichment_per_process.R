
# Load required libraries
pacman::p_load("tibble", "RColorBrewer", "pheatmap", "readr",
               "stringr","purrr", "viridis", "rvest",
               "ComplexHeatmap", "circlize", "glue","patchwork",
               "dplyr", "ggplot2", "textshape", "paletteer", "tibble", "tidyr")


# Define path to HOMER output directory
motif_dir <- "../results/RNA_ATAC_integration/Homer/pathway_enhancers_enrichment"

# List all pathway-specific motif result folders
folders <- list.dirs(motif_dir, full.names = FALSE, recursive = FALSE)

# Loop through each folder, clean the knownResults.txt, and save as CSV
for (folder in folders) {
  
  # Path to knownResults.txt
  motif_file <- glue("{motif_dir}/{folder}/knownResults.txt")
  
  # Clean output filename
  new_name <- gsub("_homer", "", folder)
  
  # Read and clean motif table
  motif_table_clean <- vroom::vroom(motif_file) %>% 
    select(`Motif Name`, `P-value`, `% of Target Sequences with Motif`) %>% 
    rename(TF_motif = `Motif Name`,
      P_value = `P-value`,
      Percent_Targets = `% of Target Sequences with Motif`) %>%
    mutate(log10_pvalue = -log10(as.numeric(P_value)),
      P_value = as.numeric(P_value),
      Percent_Targets = as.numeric(gsub("%", "", Percent_Targets)),
      TF_motif = str_extract(TF_motif, "^[^(]+"))  # Remove trailing details
  
  # Write cleaned table to CSV in same folder
  write.csv(motif_table_clean, file = glue("{motif_dir}/{folder}/{new_name}.csv"), row.names = FALSE)
}

### Load and filter cleaned motif tables for pathway-linked enhancers ----

# Define path to HOMER enrichment folder
motif_dir <- "../results/RNA_ATAC_integration/Homer/pathway_enhancers_enrichment"

# List all folders and clean names
folders <- list.dirs(motif_dir, full.names = FALSE, recursive = FALSE)
samples <- gsub("_homer", "", folders)

# Reorder samples manually for heatmap plotting
samples <- samples[c(13, 14, 1, 9, 5, 4, 6, 10, 3, 11, 2, 7, 8, 15, 12)]

# Set p-value threshold
threshold <- -log10(0.05)

# Read and filter each sample's motif enrichment table
ec_motif_list_full <- map(samples, function(sample) {
  read_csv(glue("{motif_dir}/{sample}_homer/{sample}.csv")) %>%
    mutate(pathway = sample) %>% 
    filter(log10_pvalue > threshold) %>%
    select(TF_motif, log10_pvalue, pathway) %>%
    distinct(TF_motif, .keep_all = TRUE)
})

names(ec_motif_list_full) <- samples


ec_motif_full <- bind_rows(ec_motif_list_full)


ec_motif_full_f <- ec_motif_full %>% 
  filter(pathway %in% c("Transcription_factors", "Transporters", "Angiogenesis", "Autophagy", 
                        "Hypoxia", "Coagulation", "chromatin_remodeling", "Inflammatory_response", "EC_metabolism")) %>% 
  mutate(pathway = str_replace_all(pathway, c("EC_metabolism" = "Cell metabolism", "BBB" = "Blood brain barrier", "_" = " ", "chromatin" = "Chromatin")))

ec_motif_full_wide <- ec_motif_full_f %>%  
  pivot_wider( id_cols = TF_motif,                
               names_from = pathway,              
               values_from = log10_pvalue, values_fill = 0) %>% 
  column_to_rownames("TF_motif")


ec_motif_all_pathways <-  ec_motif_full_wide %>%
  filter(rowSums(across(everything(), ~ . > 0)) >= 7)

ec_motif_all_pathways_filt <- ec_motif_all_pathways %>%
  mutate(avg_value = rowMeans(across(everything()), na.rm = TRUE),
         num_process = rowSums(across(everything(), ~ .x > 0))) %>% 
  rownames_to_column("TF_motif") %>% 
  filter(!TF_motif %in% c("EWS:ERG-fusion", "EWS:FLI1-fusion", "Unknown-ESC-element", "Ets1-distal")) %>% 
  mutate(split = case_when(str_detect(TF_motif, regex("sox", ignore_case = TRUE)) ~ "Sox",
                           str_detect(TF_motif, regex("jun|fra|fos|atf|maf|batf|ap-1|nfe2", ignore_case = TRUE)) ~ "AP-1",
                           str_detect(TF_motif, regex("etv|elk|spib|spdef|ets|pu|erg|elf|gabp|fli|ehf", ignore_case = TRUE)) ~ "Ets",
                           str_detect(TF_motif, regex("klf1$|klf3$|klf4$", ignore_case = TRUE)) ~ "Klf",
                           str_detect(TF_motif, regex("fox", ignore_case = TRUE)) ~ "Fox",
                           str_detect(TF_motif, regex("tead", ignore_case = TRUE)) ~ "Tead",
                           str_detect(TF_motif, regex("mef2", ignore_case = TRUE)) ~ "Mef2",
                           str_detect(TF_motif, regex("nfkb", ignore_case = TRUE)) ~ "Nfkb",
                           str_detect(TF_motif, regex("lef|tcf", ignore_case = TRUE)) ~ "Tcf",
                           str_detect(TF_motif, regex("gata", ignore_case = TRUE)) ~ "Gata",
                           TRUE ~ "remove")) %>% 
  filter(!split == "remove") %>% 
  group_by(split) %>% 
  arrange(desc(num_process)) %>% 
  slice_head(n = 3) %>% 
  select(-avg_value, -num_process) %>% 
  mutate(split = factor(split, level = c("AP-1", "Ets", "Sox", "Nfkb", "Klf", "Gata", "Tcf", "Fox", "Tead"))) %>% 
  column_to_rownames("TF_motif")

col_dend = dendsort::dendsort(hclust(dist(t(ec_motif_all_pathways))),  isReverse = TRUE)
#row_dend = dendsort::dendsort(hclust(dist(ec_motif_all_pathways_filt[ ,-c(10,11)])),  isReverse = TRUE)

hmap_motfs_all_pathways <- ComplexHeatmap::pheatmap(mat = ec_motif_all_pathways_filt[-c(10,11)],
                                                    name = "-log10\nP-value",
                                                    cluster_rows = FALSE, 
                                                    treeheight_col = 15,
                                                    cluster_cols = col_dend, 
                                                    cutree_cols = 3,
                                                    split = ec_motif_all_pathways_filt$split,
                                                    show_colnames = FALSE, 
                                                    angle_col = "45", 
                                                    show_rownames = TRUE,
                                                    heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),   # Legend title appearance
                                                                                labels_gp = gpar(fontsize = 7, col = "black"),
                                                                                legend_height = unit(1.5, "cm")),       # Legend labels appearance
                                                    annotation_names_col = FALSE,
                                                    border = "black",
                                                    row_title_gp = gpar(fontsize = 9, fontface = "bold"),
                                                    fontsize_row = 7.5, 
                                                    fontsize_col = 7, main = NULL,
                                                    breaks = seq(0, max(ec_motif_full_wide), length.out = 101),
                                                    color = c("white", rocket(150, begin = 1, end = 0)[1:99]),
                                                    column_title = gt_render("<span style='font-size:12pt'>TF known motif enrichement analysis</span><br>
                                                                             <span style='font-size:10pt'>Core TFs for BEC function</span>"),
                                                    border_color = "#FFFFFF")

hmap_motfs_all_pathways


# Select top 15 motifs per cell type based on log10_pvalue

ec_motif_some_wide <- ec_motif_full_wide[!rownames(ec_motif_full_wide) %in% rownames(ec_motif_all_pathways), ]

ec_motif_some_wide_f <- ec_motif_some_wide %>% 
  rownames_to_column("TF_motif") %>% 
  filter(!TF_motif %in% c("EWS:ERG-fusion", "EWS:FLI1-fusion", "Unknown-ESC-element", "Ets1-distal")) %>% 
  mutate(split = case_when(str_detect(TF_motif, regex("sox", ignore_case = TRUE)) ~ "Sox",
                           str_detect(TF_motif, regex("jun|fra|fos|atf|maf|batf|ap-1|nfe2", ignore_case = TRUE)) ~ "AP-1",
                           str_detect(TF_motif, regex("etv|elk|spib|spdef|ets|pu|erg|elf|gabp|fli|ehf", ignore_case = TRUE)) ~ "Ets",
                           str_detect(TF_motif, regex("klf1$|klf3$|klf4$", ignore_case = TRUE)) ~ "Klf",
                           str_detect(TF_motif, regex("fox", ignore_case = TRUE)) ~ "Fox",
                           str_detect(TF_motif, regex("tead", ignore_case = TRUE)) ~ "Tead",
                           str_detect(TF_motif, regex("mef2", ignore_case = TRUE)) ~ "Mef2",
                           str_detect(TF_motif, regex("nfkb", ignore_case = TRUE)) ~ "Nfkb",
                           str_detect(TF_motif, regex("lef|tcf", ignore_case = TRUE)) ~ "Tcf",
                           str_detect(TF_motif, regex("gata", ignore_case = TRUE)) ~ "Gata",
                           TRUE ~ TF_motif)) %>% 
  filter(!split %in% c("AP-1", "Ets", "Sox", "Nfkb", "Klf", "Gata", "Tcf", "Fox", "Tead")) %>% 
  select(-split) %>% 
  column_to_rownames("TF_motif") 

ec_motif_top15 <- ec_motif_some_wide_f %>%
  rownames_to_column("TF_motif") %>% 
  pivot_longer(cols = -TF_motif, values_to = "log10_pvalue", names_to = "pathway") %>% 
  group_by(TF_motif) %>% 
  arrange(desc(log10_pvalue), .by_group = T) %>% 
  distinct(TF_motif, .keep_all = T) %>% 
  ungroup() %>% 
  group_by(pathway) %>%
  slice_max(order_by = log10_pvalue, n = 4, with_ties = FALSE) %>%  # Top 10 without ties
  ungroup()

order_ec_motif_top15 <- ec_motif_top15 %>%
  group_by(TF_motif) %>%
  arrange(desc(log10_pvalue), .by_group = TRUE) %>%
  distinct(TF_motif, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(pathway = factor(pathway, levels = c("Angiogenesis", "Inflammatory response",
                                              "Coagulation", "Autophagy", "Cell metabolism",
                                              "Transporters", "Hypoxia", "Chromatin remodeling",
                                              "Transcription factors"))) %>%
  arrange(pathway, desc(log10_pvalue)) %>%
  pull("TF_motif")

row_dend = dendsort::dendsort(hclust(dist(ec_motif_some_wide_f[order_ec_motif_top15, ])),  isReverse = TRUE)
hmap_motfs_some_pathways <- ComplexHeatmap::pheatmap(mat = ec_motif_some_wide_f[order_ec_motif_top15, ],
                                                     name = "-log10\nP-value",
                                                     cluster_rows = row_dend,
                                                     treeheight_col = 15,
                                                     show_column_dend = FALSE,
                                                     clustering_method = "ward.D2",
                                                     cluster_cols = col_dend, 
                                                     cutree_cols = 3, 
                                                     show_row_dend = FALSE,
                                                     show_colnames = TRUE, main = "Predicted process-specific TF",
                                                     angle_col = "90", 
                                                     show_rownames = TRUE,
                                                     heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),   # Legend title appearance
                                                                                 labels_gp = gpar(fontsize = 7, col = "black"),
                                                                                 legend_height = unit(1, "cm")),       # Legend labels appearance
                                                    
                                                     border = "black", 
                                                     row_title = "Complementary transcription factors",
                                                     fontsize_row = 7.5, 
                                                     fontsize_col = 9, 
                                                     fontsize = 8,
                                                     row_title_gp = gpar(fontsize = 9, col = "white", fontface = "bold"),
                                                     breaks = seq(0, max(ec_motif_full_wide), length.out = 101),
                                                     color = c("white", rocket(150, begin = 1, end = 0)[1:99]),
                                                     border_color = "#FFFFFF")


png(filename = "../results/RNA_ATAC_integration/graphs/hmap_motfs_all_pathways.png", width = 3.55, height = 4, unit = "in", res = 300)
ComplexHeatmap::draw(hmap_motfs_all_pathways, show_heatmap_legend = F)
dev.off()

png(filename = "../results/RNA_ATAC_integration/graphs/hmap_motfs_complementary_pathways.png", width = 3.9, height = 6, unit = "in", res = 300)
ComplexHeatmap::draw(hmap_motfs_some_pathways, merge_legends = T, heatmap_legend_side  = "right")
dev.off()

svg(filename = "../results/RNA_ATAC_integration/graphs/hmap_motfs_combined.svg", width = 3.8, height = 8.5)
draw(hmap_motfs_all_pathways %v% hmap_motfs_some_pathways, 
     merge_legends = TRUE, gap = unit(8, "mm"), heatmap_legend_side  = "right")
dev.off()


# === Log R session info for reproducibility ===
sink("../results/logs/sessioninfo_14_heatmaps_TF_enrichment_per_process.txt")
sessionInfo()
sink()
