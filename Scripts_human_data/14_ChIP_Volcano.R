
pacman::p_load("tibble", "RColorBrewer", "pheatmap", "circlize", "stringr", "cowplot", #"gt",
               "dplyr", "ggplot2", "textshape", "paletteer", "EnhancedVolcano", "cowplot")

#### JUNB volcanos----

files <- list.files(path = "../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/deg_tables/", recursive = F)
files 

samples <- gsub(".csv", "", files)
samples

# graph folder
if (!dir.exists("../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/graphs_ChIP_JUNB/volcano")){
  dir.create("../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/graphs_ChIP_JUNB/volcano")
}

#volcano plots----

for (loop_var in samples) {
  data <- read.csv(file = sprintf("../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/deg_tables/%s.csv",loop_var))
  data <- data %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var)))) %>% 
    mutate(diffexpressed = ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) >= 1, "Up", 
                                  ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) <= -1, "Down", "No"))) %>% 
    dplyr::select(peaks, starts_with("logFC"), starts_with("FDR"), diffexpressed) %>% 
    column_to_rownames("peaks")
  
  assign(paste0(loop_var), value = data)
}


full_names <- gsub("CCM", "(CCM env)", gsub("VEH", "(Veh)", gsub("_", " ", samples)))


names(full_names) <- samples
full_names

for (loop_var in samples) {
  data <- get(loop_var)
  summary_data <- data %>% 
    group_by(diffexpressed) %>% 
    summarise(n = n())
  
  full_name <- full_names[loop_var]
  
  title <- paste0("ChIP JUNB - DAR\n", full_name)
  
  volcano_p <- EnhancedVolcano(data, 
                               lab = NA,
                               legendLabSize = 12,
                               axisLabSize = 13,
                               x = sprintf("logFC_%s", loop_var),
                               y = sprintf("FDR_%s", loop_var), 
                               title = title, subtitleLabSize = 17,
                               subtitle = paste0(summary_data[1,]$n, " lost and ", summary_data[3,]$n, " gained"),
                               pCutoff = 0.05,
                               FCcutoff = 1, 
                               raster = T) 
  ggsave(filename = paste0("../results/Human_cells/ChipSeq/edgeR_results_ChIP_JUNB/graphs_ChIP_JUNB/volcano/", loop_var, "_volcano.png"),
         plot = volcano_p, width = 6, height = 6, units = "in", dpi = 300, device = "png")
  #assign(paste0(loop_var, "_hmap"), value = hm)
  
}

#### H3K27ac volcanos----

files <- list.files(path = "../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/deg_tables/", recursive = F)
files 

samples <- gsub(".csv", "", files)
samples

# graph folder
if (!dir.exists("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/graphs_ChIP_H3K27ac/volcano")){
  dir.create("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/graphs_ChIP_H3K27ac/volcano")
}

#volcano plots----

for (loop_var in samples) {
  data <- read.csv(file = sprintf("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/deg_tables/%s.csv",loop_var))
  data <- data %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var)))) %>% 
    mutate(diffexpressed = ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) >= 1, "Up", 
                                  ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) <= -1, "Down", "No"))) %>% 
    dplyr::select(peaks, starts_with("logFC"), starts_with("FDR"), diffexpressed) %>% 
    column_to_rownames("peaks")
  
  assign(paste0(loop_var), value = data)
}


full_names <- gsub("CCM", "(CCM env)", gsub("VEH", "(Veh)", gsub("_", " ", samples)))


names(full_names) <- samples
full_names

for (loop_var in samples) {
  data <- get(loop_var)
  summary_data <- data %>% 
    group_by(diffexpressed) %>% 
    summarise(n = n())
  
  full_name <- full_names[loop_var]
  
  title <- paste0("ChIP H3K27ac - DAR\n", full_name)
  
  volcano_p <- EnhancedVolcano(data, 
                               lab = NA,
                               legendLabSize = 12,
                               axisLabSize = 13,
                               x = sprintf("logFC_%s", loop_var),
                               y = sprintf("FDR_%s", loop_var), 
                               title = title, subtitleLabSize = 17,
                               subtitle = paste0(summary_data[1,]$n, " lost and ", summary_data[3,]$n, " gained"),
                               pCutoff = 0.05,
                               FCcutoff = 1, 
                               raster = T) 
  ggsave(filename = paste0("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K27ac/graphs_ChIP_H3K27ac/volcano/", loop_var, "_volcano.png"),
         plot = volcano_p, width = 6, height = 6, units = "in", dpi = 300, device = "png")
  #assign(paste0(loop_var, "_hmap"), value = hm)
  
}

#### H3K4me1 volcanos----

files <- list.files(path = "../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/deg_tables/", recursive = F)
files 

samples <- gsub(".csv", "", files)
samples

# graph folder
if (!dir.exists("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/graphs_ChIP_H3K4me1/volcano")){
  dir.create("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/graphs_ChIP_H3K4me1/volcano")
}

#volcano plots----

for (loop_var in samples) {
  data <- read.csv(file = sprintf("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/deg_tables/%s.csv",loop_var))
  data <- data %>% 
    arrange(desc(get(sprintf("logFC_%s", loop_var)))) %>% 
    mutate(diffexpressed = ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) >= 0.5, "Up", 
                                  ifelse(get(sprintf("FDR_%s", loop_var)) <= 0.05 & get(sprintf("logFC_%s", loop_var)) <= -0.5, "Down", "No"))) %>% 
    dplyr::select(peaks, starts_with("logFC"), starts_with("FDR"), diffexpressed) %>% 
    column_to_rownames("peaks")
  
  assign(paste0(loop_var), value = data)
}


full_names <- gsub("CCM", "(CCM env)", gsub("VEH", "(Veh)", gsub("_", " ", samples)))


names(full_names) <- samples
full_names

for (loop_var in samples) {
  data <- get(loop_var)
  summary_data <- data %>% 
    group_by(diffexpressed) %>% 
    summarise(n = n()) %>% 
    complete(diffexpressed = c("Up", "Down", "No"), fill = list(n = 0))
  
  subtitle_volcano <- paste0(summary_data %>% filter(diffexpressed == "Down") %>% pull(n), " lost and ", 
    summary_data %>% filter(diffexpressed == "Up") %>% pull(n), " gained")
  
  full_name <- full_names[loop_var]
  
  title <- paste0("ChIP H3K4me1 - DAR\n", full_name)
  
  volcano_p <- EnhancedVolcano(data, 
                               lab = NA,
                               legendLabSize = 12,
                               axisLabSize = 13,
                               x = sprintf("logFC_%s", loop_var),
                               y = sprintf("FDR_%s", loop_var), 
                               title = title, subtitleLabSize = 17,
                               subtitle = paste0(summary_data[1,]$n, " lost and ", summary_data[3,]$n, " gained"),
                               pCutoff = 0.05,
                               FCcutoff = 1, 
                               raster = T) 
  ggsave(filename = paste0("../results/Human_cells/ChipSeq/edgeR_results_ChIP_H3K4me1/graphs_ChIP_H3K4me1/volcano/", loop_var, "_volcano.png"),
         plot = volcano_p, width = 6, height = 6, units = "in", dpi = 300, device = "png")
  #assign(paste0(loop_var, "_hmap"), value = hm)
  
}

# Optionally save session info
sink("../logs/sessioninfo_human_14_ChIP_volcano.txt")
sessionInfo()
sink()




