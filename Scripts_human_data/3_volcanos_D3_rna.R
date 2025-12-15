
pacman::p_load("tibble", "RColorBrewer", "pheatmap", "circlize", "stringr", "cowplot",
               "dplyr", "ggplot2", "textshape", "paletteer", "EnhancedVolcano", "cowplot")


files <- list.files(path = "edgeR_results_bulk_RNA/deg_tables/", recursive = F)
files 

samples <- gsub(".csv", "", files[-1])

# graph folder
if (!dir.exists("graphs_RNA")){
  dir.create("graphs_RNA")
}

#volcano plots----

data_volcano <- list()
for (loop_var in samples) {
  data <- read.csv(file = sprintf("edgeR_results_bulk_RNA/deg_tables/%s.csv",loop_var))
  data <- data %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var)))) %>% 
    mutate(diffexpressed = ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) >= 0.5, "Up", 
                                  ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) <= -0.5, "Down", "No"))) %>% 
    dplyr::select(gene, starts_with("logFC"), starts_with("FDR"), diffexpressed) %>% 
    column_to_rownames("gene")
  
  data_volcano[[loop_var]] <- data
}

full_names <- gsub("CCT", "(CCM env + T5224)", gsub("T52", "(T5224)", gsub("CCM", "(CCM env)", gsub("VEH", "(Veh)", gsub("_", " ", samples)))))


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
  
  volcano_p <- EnhancedVolcano(data, 
                               lab = NA,
                               legendLabSize = 12,
                               axisLabSize = 13, 
                               x = sprintf("logFC_%s", loop_var),
                               y = sprintf("FDR_%s", loop_var), 
                               title = title, subtitleLabSize = 17,
                               subtitle = paste0(summary_data[1,]$n, " downregulated and ", summary_data[3,]$n, " upregulated"),
                               pCutoff = 0.05,
                               FCcutoff = 0.5, 
                               raster = T) 
  ggsave(filename = paste0("graphs_RNA/volcano/", loop_var, "_volcano.png"),
         plot = volcano_p, width = 6, height = 6, units = "in", dpi = 300, device = "png")
  #assign(paste0(loop_var, "_hmap"), value = hm)
  volcano_plots[[loop_var]] <- volcano_p
  
}

saveRDS(volcano_plots[["KD_VEH_vs_WT_VEH"]], "graphs_RNA/volcano/KD_V_WT_V.rds")
saveRDS(volcano_plots[["KD_CCM_vs_WT_VEH"]], "graphs_RNA/volcano/KD_c_WT_V.rds")
saveRDS(volcano_plots[["WT_CCT_vs_WT_VEH"]], "graphs_RNA/volcano/WT_C_WT_V.rds")



