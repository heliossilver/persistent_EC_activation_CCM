
pacman::p_load("dplyr", "vroom", "stringr", "pheatmap", "purrr",
               "ggplot2", "tidyr", "tibble", "RColorBrewer")

#### 1. Clean motif enrichment output files from HOMER ----
# Define input and output paths
motif_dir <- "../results/Human_cells/Homer/DAR_results_05_enrichment/"

results_tables_homer <- "../results/Human_cells/Homer/results_tables/DAR_results_05_results"
dir.create(results_tables_homer, recursive = TRUE, showWarnings = FALSE)


# List cell-type specific motif result folders
folders <- list.dirs(motif_dir, full.names = FALSE, recursive = FALSE)


# Loop through each file and clean the knownResults.txt
for (folder in folders) {
  
  motif_file <- glue("{motif_dir}/{folder}/knownResults.txt")
  new_name <- paste0(folder, "_clean")
  
  motif_df <- vroom::vroom(motif_file) %>% 
    select(`Motif Name`, `P-value`, `% of Target Sequences with Motif`) %>% 
    dplyr::rename(TF_motif = `Motif Name`,
                  P_value = `P-value`,
                  Percent_Targets = `% of Target Sequences with Motif`) %>% 
    mutate(log10_pvalue = -log10(as.numeric(P_value)),
           Percent_Targets = as.numeric(gsub("%", "", Percent_Targets)),
           TF_motif = str_extract(TF_motif, "^[^(]+"))  # Remove trailing details
  
  # Write cleaned CSV
  write.csv(motif_df, file = glue("{motif_dir}/{folder}/{new_name}.csv"), row.names = FALSE)
}


# Precalculate -log10(0.05) threshold
log10_thresh <- -log10(0.05)

# Read and filter cleaned motif tables
motif_list <- map(folders, function(folder) {
  read_csv(glue("{motif_dir}/{folder}/{folder}_clean.csv")) %>%
    filter(log10_pvalue >= log10_thresh) %>%
    mutate(Classification = gsub("_homer05_up", "", folder)) %>%
    distinct(TF_motif, .keep_all = TRUE) %>%
    select(TF_motif, log10_pvalue, Classification)
})

names(motif_list) <- gsub("_homer05_up", "", folders)

ec_motif_up <- bind_rows(motif_list[c(2, 4, 6)])
ec_motif_up_enhancers <- bind_rows(motif_list[c(1, 3, 5)])

write.csv(ec_motif_up, file = glue("{results_tables_homer}/ec_motif_up_full.csv"), row.names = FALSE)

write.csv(ec_motif_up_enhancers, file = glue("{results_tables_homer}/ec_motif_up_enhancers.csv"), row.names = FALSE)

















files <- unique(gsub("_homer05_up_full|_homer05_up_enhancers", "", folders))
files <- files[c(2, 1, 3)]

### Using top 10 per cell type----

motifs_up <- c()
ec_motif_list_up <- list()
for (loop_var in files) {
  #Read data
  data.motif_up.plot <- vroom(file = sprintf("Homer/DAR_1_TF_enrichment/%s_homer01_up/knownResults.txt", loop_var)) %>% 
    mutate(TF_motif = str_extract(`Motif Name`, "^[^(]+"),
           'condition' = loop_var,
           '-log10 P-value' = -log10(`P-value`))
  
  
  data.motif_up.plot <- data.motif_up.plot %>% 
    distinct(TF_motif, .keep_all = T) %>% 
    select(TF_motif, `-log10 P-value`, condition)
  #Get the top 10 gene names
  tf_up <-  data.motif_up.plot %>% 
    head(n = 15) %>% 
    pull(TF_motif)
  
  #Put data for plotting in a list
  ec_motif_list_up[[loop_var]] <- data.motif_up.plot
  
  #get all the top gene names 
  motifs_up <- unique(c(motifs_up, tf_up))
  
  
}

motifs_down <- c()
ec_motif_list_down <- list()
for (loop_var in files) {
  #Read data
  data.motif_down.plot <- vroom(file = sprintf("Homer/DAR_1_TF_enrichment/%s_homer01_down/knownResults.txt", loop_var)) %>% 
    mutate(TF_motif = str_extract(`Motif Name`, "^[^(]+"),
           'condition' = loop_var,
           '-log10 P-value' = -log10(`P-value`))
  
  
  data.motif_down.plot <- data.motif_down.plot %>% 
    distinct(TF_motif, .keep_all = T) %>% 
    select(TF_motif, `-log10 P-value`, condition)
  
  #Get the top 10 gene names
  tf_down <-  data.motif_down.plot %>% 
    head(n = 15) %>% 
    pull(TF_motif)
  
  #Put data for plotting in a list
  ec_motif_list_down[[loop_var]] <- data.motif_down.plot
  
  #get all the top gene names 
  motifs_down <- unique(c(motifs_down, tf_down))
  
  
}


#new_levels <- gsub("_vs_Flox", "", files)


##up data for heatmap ----
ec_motif_up <- bind_rows(ec_motif_list_up)
ec_motif_up_filt <- ec_motif_up %>% 
  filter(TF_motif %in% motifs_up) 
ec_motif_up_filt$TF_motif <- factor(x = ec_motif_up_filt$TF_motif, levels = motifs_up)
ec_motif_up_filt$condition <- factor(x = ec_motif_up_filt$condition, levels = files)

ec_motif_up_wide <- ec_motif_up_filt %>%                      
  pivot_wider( id_cols = TF_motif,                
               names_from = condition,              
               values_from = `-log10 P-value`) %>% 
  column_to_rownames("TF_motif")



##down data for heatmap ----


ec_motif_down <- bind_rows(ec_motif_list_down)
ec_motif_down_filt <- ec_motif_down %>% 
  filter(TF_motif %in% motifs_down)

ec_motif_down_filt$condition <- factor(x = ec_motif_down_filt$condition, levels = files)



ec_motif_down_wide <- ec_motif_down_filt %>%  
  pivot_wider( id_cols = TF_motif,                
               names_from = condition,              
               values_from = `-log10 P-value`) %>% 
  column_to_rownames("TF_motif")


#annotations for heatmaps and heatmaps ----

annot <- data.frame(Condition = str_replace_all(files, c("CCM" = "(CCM env)", "T52" = "(T5224)",
                                                         "CCT" = "(CCM env + T5224)", "ATAC_" = "",
                                                         "_" = " ", "VEH" = "(Veh)", " vs " = "\nvs ")),
                    row.names = colnames(ec_motif_up_wide))
colnames(annot) <- "Condition"

annot$Condition <- factor(annot$Condition, levels = unique(annot$Condition))



# Generate the color palette
cond_colors <- c("#FB6A4A", "#756BB1" , "#3182BD", "#31A354")
names(cond_colors) <-  str_replace_all(files, c("CCM" = "(CCM env)", "T52" = "(T5224)",
                                                "CCT" = "(CCM env + T5224)", "ATAC_" = "",
                                                "_" = " ", "VEH" = "(Veh)", " vs " = "\nvs "))

my_colour <-  list("Condition" = cond_colors)

# ## Heatmaps 


ec_motif_up_wide[ec_motif_up_wide == -0] <- 0
hmap_up05 <- ComplexHeatmap::pheatmap( mat = ec_motif_up_wide, 
                                       name = "-log10 P-Value", 
                                       cluster_rows = F, cluster_cols = F,
                                       show_colnames = F, fontsize_col = 7, 
                                       breaks = seq(0, 400, length.out = 100),
                                       color = colorRampPalette(brewer.pal(9, "Reds"))(100),
                                       main = "TF motif enrichment from open\ncCRE (Positive cCRE-RNA corr.)",
                                       annotation_col = annot, #angle_col = "45",
                                       display_numbers = T,number_format = "%1.0f", number_color = "black",
                                       annotation_colors = my_colour)    

hmap_up05_cl <- ComplexHeatmap::pheatmap( mat = ec_motif_up_wide, 
                                          name = "-log10 P-Value", 
                                          cluster_rows = T, cluster_cols = F, treeheight_row = 25,
                                          show_colnames = F, fontsize_col = 7,
                                          breaks = seq(0, 400, length.out = 100),
                                          color = colorRampPalette(brewer.pal(9, "Reds"))(100),
                                          main = "TF motif enrichment from open\ncCRE (Positive cCRE-RNA corr.)",
                                          annotation_col = annot, #angle_col = "45",
                                          display_numbers = T,number_format = "%1.0f", number_color = "black",
                                          annotation_colors = my_colour)  



ec_motif_down_wide[ec_motif_down_wide == -0] <- 0
hmap_down05 <- ComplexHeatmap::pheatmap(mat = ec_motif_down_wide, 
                                        name = "-log10 P-Value", 
                                        cluster_rows = F, cluster_cols = F, 
                                        show_colnames = F, fontsize_col = 7,
                                        display_numbers = T,number_format = "%1.0f", number_color = "black",
                                        breaks = seq(0, 400, length.out = 100),
                                        color = colorRampPalette(brewer.pal(9, "Blues"))(100),
                                        main = "Known motif enrichment from\n lost CA (logFC<1 & FDR<0.01)",
                                        annotation_col = annot, 
                                        annotation_colors = my_colour)    

hmap_down05_cl <- ComplexHeatmap::pheatmap( mat = ec_motif_down_wide, 
                                            name = "-log10 P-Value", display_numbers = T, 
                                            cluster_rows = T, cluster_cols = F, treeheight_row = 25,
                                            show_colnames = F, fontsize_col = 7,
                                            number_format = "%1.0f", number_color = "black",
                                            breaks = seq(0, 400, length.out = 100),
                                            color = colorRampPalette(brewer.pal(9, "Blues"))(100),
                                            main = "Known motif enrichment from\n lost CA (logFC<1 & FDR<0.01)",
                                            annotation_col = annot, #angle_col = "45", 
                                            annotation_colors = my_colour) 





if (!dir.exists("Homer/homer_graphs/FDR_01")) {
  dir.create("Homer/homer_graphs/FDR_01", recursive = T)
}

png(filename = "Homer/homer_graphs/FDR_01/dar01up_KO_top20.png", width = 5.5, height = 5, units = "in", res = 300)
ComplexHeatmap::draw(hmap_up05, merge_legends = T)
dev.off()

png(filename = "Homer/homer_graphs/FDR_01/dar01down_KO_top20.png", width = 5.8, height = 5.5, units = "in", res = 300)
ComplexHeatmap::draw(hmap_down05, merge_legends = T)
dev.off()


png(filename = "Homer/homer_graphs/FDR_01/dar01up_KO_top20_cluster.png", width = 6, height = 5.5, units = "in", res = 300)
ComplexHeatmap::draw(hmap_up05_cl, merge_legends = T)
dev.off()

png(filename = "Homer/homer_graphs/FDR_01/dar01down_KO_top20_cluster.png", width = 6, height = 5.5, units = "in", res = 300)
ComplexHeatmap::draw(hmap_down05_cl, merge_legends = T)
dev.off()



hmap_up_down_list <- list(grid::grid.grabExpr(ComplexHeatmap::draw(hmap_up05, show_heatmap_legend = T, show_annotation_legend = F)),
                          grid::grid.grabExpr(ComplexHeatmap::draw(hmap_down05, merge_legends = T)))

hmap_up_down_plots <- patchwork::wrap_plots(hmap_up_down_list, ncol = 2, widths = c(1, 1.1))
ggsave(filename = "Homer/homer_graphs/FDR_01/hmap_TF_updown_FDR_01.png", plot = hmap_up_down_plots, device = "png", width = 9.2, height = 6, units = "in", dpi = 300)


hmap_up_down_cluster_list <- list(grid::grid.grabExpr(ComplexHeatmap::draw(hmap_up05_cl, show_heatmap_legend = T, show_annotation_legend = FALSE)),
                                  grid::grid.grabExpr(ComplexHeatmap::draw(hmap_down05_cl, merge_legends = T)))

hmap_up_down_cluster_plots <- patchwork::wrap_plots(hmap_up_down_cluster_list, ncol = 2, widths = c(1, 1.1))
ggsave(filename = "Homer/homer_graphs/FDR_01/hmap_TF_updown_clust_FDR_01.png", plot = hmap_up_down_cluster_plots, device = "png", width = 9, height = 6, units = "in", dpi = 300)

##tables ----

if (!dir.exists("Homer/results_tables")) {
  dir.create("Homer/results_tables")
}




