#install.packages("devtools")
# library(devtools)
# install_github("wjawaid/enrichR")

pacman::p_load("enrichR", "dplyr", "tibble", "tidyr", "glue", "forcats",
               "ggplot2", "vroom", "stringr", "viridis", "RColorBrewer")

# ==== SETUP ====

# Set paths and database
rna_dir <- "../results/Human_cells/D3_RNA/edgeR/deg_difexp_tables"
atac_dir <- "../results/Human_cells/D3_ATAC/edgeR/dar_difacc_tables"
background_file <- "../results/Human_cells/D3_RNA/data/log2cpm_filt.csv"
kegg_folder <- "../results/Human_cells/D3_ATAC/tables"
dir.create(kegg_folder, showWarnings = FALSE, recursive = TRUE)
irrelevant_terms <- read.csv(file = "utils/Irrelevant_terms.csv") %>% pull("Irrelevant_terms")

# Load enhancer–gene links
pcc_filt <- read.csv("../results/Human_cells/D3_ATAC/data/pcc_ccre_genes_500kb_filt.csv") %>% select(peaks, gene)


# Set Enrichr site
setEnrichrSite("Enrichr")

dbs <- listEnrichrDbs()
dbs <- "KEGG_2021_Human"


### EnrichR pipe ----
###RNA----
rna_files <- list.files(rna_dir, pattern = "_diffexp.csv", full.names = TRUE)

# Sort into WT vs Rest and KO vs Flox
samples <- rna_files[grepl("_vs_WT_VEH", rna_files)][c(4, 1, 5)]

# ==== FUNCTION TO LOAD AND FORMAT ====

read_diffexp <- function(file_path) {
  sample_name <- gsub("_vs.*$", "", basename(file_path))
  vroom(file_path) %>%
    dplyr::select(gene, diffexpressed) %>%
    mutate(sample = sample_name)
}


# Apply to all files
df_g_list <- lapply(samples, read_diffexp)
# Create named lists
names(df_g_list) <- sapply(df_g_list, function(x) unique(x$sample))


# ==== BACKGROUND GENE SET ====
bg_genes <- read.csv(background_file) %>% pull(gene)

# ==== SPLIT UP/DOWN GENES ====

get_genes_by_direction <- function(df_list, direction = "Up") {
  lapply(df_list, function(df) {
    df %>%
      filter(diffexpressed == direction) %>%
      pull(gene)
  })
}

genes_up   <- get_genes_by_direction(df_g_list, "Up")
genes_down <- get_genes_by_direction(df_g_list, "Down")

## ==== ATAC ====
# List differential expression result files
atac_files <- list.files(atac_dir, pattern = "_diffacc.csv", full.names = TRUE)

# Sort into WT vs Rest and KO vs Flox
samples_p <- atac_files[grepl("_vs_WT_VEH", atac_files)][c(4, 1, 5)]

# ==== FUNCTION TO LOAD AND FORMAT ====

read_diffacc <- function(file_path) {
  sample_name <- gsub("_vs.*$", "", basename(file_path))
  vroom(file_path) %>%
    select(peaks, diffaccessible) %>%
    left_join(pcc_filt, by = "peaks") %>%
    filter(!is.na(gene)) %>%  # Keep only linked peaks
    mutate(sample = sample_name)
  
}

# Apply to all WT and KO files
df_p_list <- lapply(samples_p, read_diffacc)
names(df_p_list) <- sapply(df_p_list, function(x) unique(x$sample))


# ==== BACKGROUND GENE SET for linked genes ====
bg_atac <- unique(pcc_filt$gene)

# ==== SPLIT UP/DOWN GENES FROM PEAK DATA ====

get_genes_peaks_by_direction <- function(df_list, direction = "Up") {
  lapply(df_list, function(df) {
    df %>%
      filter(diffaccessible == direction) %>%
      pull(gene)
  })
}

peaks_up   <- get_genes_peaks_by_direction(df_p_list, "Up")
peaks_down <- get_genes_peaks_by_direction(df_p_list, "Down")

# ==== RUN ENRICHR ====

run_enrichr_list <- function(gene_lists, db, background = NULL) {
  lapply(gene_lists, function(gene_set) {
    enrichr(gene_set, db, background, include_overlap = TRUE)
  })
}

enrichr_g_up   <- run_enrichr_list(genes_up, dbs, bg_genes)
enrichr_g_down <- run_enrichr_list(genes_down, dbs, bg_genes)

enrichr_p_up   <- run_enrichr_list(peaks_up, dbs, bg_atac)
enrichr_p_down <- run_enrichr_list(peaks_down, dbs, bg_atac)


# Function to process Enrichr results
process_file <- function(enrichr_list, condition, status, source) {
  if (!condition %in% names(enrichr_list)) return(NULL)
  
  enrich_df <- as.data.frame(enrichr_list[[condition]]$KEGG_2021_Human)
  
  enrich_df %>% 
    mutate(status = status, 
           source = source,
           condition = condition) %>% 
    filter(P.value <= 0.05,
           !Term %in% irrelevant_terms) %>% 
    separate(col = Overlap, into = c("Count", "total"), sep = "/", convert = T, remove = FALSE) %>% 
    mutate("log10_AdjPval" = -log10(Adjusted.P.value),
           "log10_Pval" = -log10(P.value)) %>%
    select(-contains("Old"))
  
}


conditions <-  Reduce(intersect, list(names(enrichr_g_up), names(enrichr_g_down),
                                         names(enrichr_p_up), names(enrichr_p_down)))

enrichR_list <- list()
for (loop_var in conditions) {
  
  #bind data
  data_wt <- bind_rows(process_file(enrichr_g_up, loop_var, status = "Activated",source =  "Genes"),
                       process_file(enrichr_g_down, loop_var, "Suppressed", "Genes"),
                       process_file(enrichr_p_up, loop_var, "Activated", "Enhancers"),
                       process_file(enrichr_p_down, loop_var, "Suppressed", "Enhancers"))
  
  enrichR_list[[loop_var]] <- data_wt
  
}

# Combine all WT and KO results into long data frames
enrichR_df <- bind_rows(enrichR_list, .id = "condition") %>% 
  mutate(condition = factor(condition, levels = conditions),
         source = factor(source, levels = c("Genes", "Enhancers")))

write.csv(enrichR_df, "../results/Human_cells/D3_ATAC/tables/KEGG_enrichment.csv", row.names = FALSE)


top_terms <- enrichR_df %>%
  filter(status == "Activated") %>% 
  group_by(condition, source) %>%
  slice_max(order_by = log10_Pval, n = 10) %>% 
  pull(Term) %>% 
  unique()

new_names <- c("siPDCD10", "siPDCD10 + CCM-like env", "Control + CCM-like env")
names(new_names) <- levels(enrichR_df$condition)

data_plot <- enrichR_df %>% 
  mutate(condition = recode(as.character(condition), !!!new_names)) %>% 
  filter(Term %in% top_terms,
         status == "Activated") %>% 
  mutate(condition = factor(condition, levels = new_names)) %>% 
  arrange(source, condition, desc(log10_Pval)) %>% 
  mutate(Term = fct_rev(fct_inorder(Term))) 


cond_colors <- c("#FB9A99", "#E31A1C", "#1F78B4")


plot_dot_kegg <- function(data) {
  ggplot(data, aes(x = condition, y = Term, color = log10_Pval, size = Count)) +
    geom_point() +
    scale_color_viridis_c(breaks = seq(floor(min(data$log10_Pval)), ceiling(max(data$log10_Pval)), length.out = 5),
                          labels = function(x) round(x, 0),name = "-log10\nP-value") +
    theme_bw() +
    ggtitle("KEGG enrichment") +
    theme(axis.text.y = element_text(size = 9, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "bold", colour = cond_colors),
          legend.key.height = unit(12, "pt"),
          legend.title = element_text(size = 10, face = "bold", margin = margin(t = 2, b = 4)),
          legend.text = element_text(size = 9),
          axis.ticks.x = element_line(linewidth = 0),
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 9.5, face = "bold", hjust = 0.5, vjust = 0.5, colour = "black")) +
    scale_size_continuous(range = c(1, 3)) +
    facet_wrap(~source, scales = "free_x")
}


plot_dot_kegg(data_plot)


selected_terms <- read.csv("utils/selected_terms.csv") %>% pull("selected_pathways")
selected_terms <- selected_terms[-c(11, 19, 12, 20)]
cond_colors <- c( "#E31A1C", "#1F78B4")

selected_to_plot <- enrichR_df %>% 
  filter(condition %in% c("KD_CCM", "WT_CCM")) %>% 
  mutate(condition = recode(as.character(condition), !!!new_names)) %>% 
  filter(Term %in% selected_terms,
         status == "Activated") %>% 
  mutate(condition = factor(condition, levels = new_names)) %>% 
  arrange(source, condition, desc(log10_Pval)) %>% 
  mutate(Term = fct_rev(fct_inorder(Term))) 


plot_dot_kegg(selected_to_plot)

ggsave("../results/Human_cells/plots/kegg_Selected.png", plot_dot_kegg(selected_to_plot), height = 3.8, width = 4.7, units = "in", dpi = 300)


# Optionally save session info
sink("../logs/sessioninfo_human_12_Enrichr_KEGG.txt")
sessionInfo()
sink()
