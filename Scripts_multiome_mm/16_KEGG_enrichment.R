# Load required packages
pacman::p_load("enrichR", "dplyr", "tibble", "tidyr", "glue", "forcats", 
               "ggplot2", "vroom", "stringr", "viridis", "RColorBrewer")

# ==== SETUP ====

# Set paths and database
rna_dir <- "../results/RNA/edgeR/deg_difexp_tables/"
atac_dir <- "../results/ATAC/edgeR/dar_diffacc_tables/"
background_file <- "../results/RNA/data/log2cpm_filtered.csv"
irrelevant_terms <- read.csv(file = "utils/Irrelevant_terms.csv") %>% pull("Irrelevant_terms")

# Load enhancer–gene links
links_df <- read.csv("../results/ATAC/data/endo_linkpeaks.csv")
links_to_join <- links_df %>%
  select(gene, peaks)

# Set Enrichr site
setEnrichrSite("Enrichr")
dbs <- "KEGG_2019_Mouse"

# ==== LOAD FILES ====

## ==== RNA ====
# List differential expression result files
rna_files <- list.files(rna_dir, pattern = "_diffexp.csv", full.names = TRUE)

# Sort into WT vs Rest and KO vs Flox
samples_wt_g <- rna_files[grepl("_vs_rest", rna_files)][c(1, 3, 2, 4, 5)]
samples_ko_g <- rna_files[grepl("_vs_Flox", rna_files)][c(1, 3, 2, 4, 5)]

# ==== FUNCTION TO LOAD AND FORMAT ====

read_diffexp <- function(file_path) {
  sample_name <- gsub("_vs.*$", "", basename(file_path))
  vroom(file_path) %>%
    select(gene, diffexpressed) %>%
    mutate(sample = sample_name)
}

# Apply to all WT and KO files
df_wt_g_list <- lapply(samples_wt_g, read_diffexp)
df_ko_g_list <- lapply(samples_ko_g, read_diffexp)

# Create named lists
names(df_wt_g_list) <- sapply(df_wt_g_list, function(x) unique(x$sample))
names(df_ko_g_list) <- sapply(df_ko_g_list, function(x) unique(x$sample))

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

genes_wt_up   <- get_genes_by_direction(df_wt_g_list, "Up")
genes_wt_down <- get_genes_by_direction(df_wt_g_list, "Down")

genes_ko_up   <- get_genes_by_direction(df_ko_g_list, "Up")
genes_ko_down <- get_genes_by_direction(df_ko_g_list, "Down")


## ==== ATAC ====
# List differential expression result files
atac_files <- list.files(atac_dir, pattern = "_diffacc.csv", full.names = TRUE)

# Sort into WT vs Rest and KO vs Flox
samples_wt_p <- atac_files[grepl("_vs_rest", atac_files)][c(1, 3, 2, 4, 5)]
samples_ko_p <- atac_files[grepl("_vs_Flox", atac_files)][c(1, 3, 2, 4, 5)]

# ==== FUNCTION TO LOAD AND FORMAT ====

read_diffacc <- function(file_path) {
  sample_name <- gsub("_vs.*$", "", basename(file_path))
  vroom(file_path) %>%
    select(peaks, diffaccessible) %>%
    left_join(links_to_join, by = "peaks") %>%
    filter(!is.na(gene)) %>%  # Keep only linked peaks
    mutate(sample = sample_name)
  
}

# Apply to all WT and KO files
df_wt_p_list <- lapply(samples_wt_p, read_diffacc)
df_ko_p_list <- lapply(samples_ko_p, read_diffacc)

# Create named lists
names(df_wt_p_list) <- sapply(df_wt_p_list, function(x) unique(x$sample))
names(df_ko_p_list) <- sapply(df_ko_p_list, function(x) unique(x$sample))

# ==== BACKGROUND GENE SET for linked genes ====
bg_atac <- unique(links_to_join$gene)
# ==== SPLIT UP/DOWN GENES FROM PEAK DATA ====

get_genes_peaks_by_direction <- function(df_list, direction = "Up") {
  lapply(df_list, function(df) {
    df %>%
      filter(diffaccessible == direction) %>%
      pull(gene)
  })
}

peaks_wt_up   <- get_genes_peaks_by_direction(df_wt_p_list, "Up")
peaks_wt_down <- get_genes_peaks_by_direction(df_wt_p_list, "Down")

peaks_ko_up   <- get_genes_peaks_by_direction(df_ko_p_list, "Up")
peaks_ko_down <- get_genes_peaks_by_direction(df_ko_p_list, "Down")

# ==== RUN ENRICHR ====

run_enrichr_list <- function(gene_lists, db, background) {
  lapply(gene_lists, function(gene_set) {
    enrichr(gene_set, db, background, include_overlap = TRUE)
  })
}

enrichr_wt_g_up   <- run_enrichr_list(genes_wt_up, dbs, bg_genes)
enrichr_wt_g_down <- run_enrichr_list(genes_wt_down, dbs, bg_genes)

enrichr_ko_g_up   <- run_enrichr_list(genes_ko_up, dbs, bg_genes)
enrichr_ko_g_down <- run_enrichr_list(genes_ko_down, dbs, bg_genes)

enrichr_wt_p_up   <- run_enrichr_list(peaks_wt_up, dbs, bg_atac)
enrichr_wt_p_down <- run_enrichr_list(peaks_wt_down, dbs, bg_atac)

enrichr_ko_p_up   <- run_enrichr_list(peaks_ko_up, dbs, bg_atac)
enrichr_ko_p_down <- run_enrichr_list(peaks_ko_down, dbs, bg_atac)


# Function to process Enrichr results
process_file <- function(enrichr_list, cell_type, status, source) {
  if (!cell_type %in% names(enrichr_list)) return(NULL)
  
  enrich_df <- as.data.frame(enrichr_list[[cell_type]]$KEGG_2019_Mouse)
  
  enrich_df %>% 
    mutate(status = status, 
           source = source,
           cell_type = cell_type) %>% 
    filter(P.value <= 0.05,
           !Term %in% irrelevant_terms) %>% 
    separate(col = Overlap, into = c("Count", "total"), sep = "/", convert = T, remove = FALSE) %>% 
    mutate("log10_AdjPval" = -log10(Adjusted.P.value),
           "log10_Pval" = -log10(P.value)) %>%
    select(-contains("Old"))
  
}

cell_types_wt <-  Reduce(intersect, list(names(enrichr_wt_g_up), names(enrichr_wt_g_down),
                                   names(enrichr_wt_p_up), names(enrichr_wt_p_down)))
cell_types_ko <- Reduce(intersect, list(names(enrichr_ko_g_up), names(enrichr_ko_g_down),
                                 names(enrichr_ko_p_up), names(enrichr_ko_p_down)))

enrichR_list_wt <- list()
for (loop_var in cell_types_wt) {
  
  #bind data
  data_wt <- bind_rows(process_file(enrichr_wt_g_up, loop_var, status = "Activated",source =  "Genes"),
                       process_file(enrichr_wt_g_down, loop_var, "Suppressed", "Genes"),
                       process_file(enrichr_wt_p_up, loop_var, "Activated", "Enhancers"),
                       process_file(enrichr_wt_p_down, loop_var, "Suppressed", "Enhancers"))
  
  enrichR_list_wt[[loop_var]] <- data_wt

}

enrichR_list_ko <- list()
for (loop_var in cell_types_ko) {
 
  data_ko <- bind_rows(process_file(enrichr_ko_g_up, loop_var, "Activated", "Genes"),
                       process_file(enrichr_ko_g_down, loop_var, "Suppressed", "Genes"),
                       process_file(enrichr_ko_p_up, loop_var, "Activated", "Enhancers"),
                       process_file(enrichr_ko_p_down, loop_var, "Suppressed", "Enhancers"))
  
  enrichR_list_ko[[loop_var]] <- data_ko
  
}


# Combine all WT and KO results into long data frames
enrichR_df_wt <- bind_rows(enrichR_list_wt, .id = "cell_type") %>% 
  mutate(cell_type = factor(cell_type, levels = cell_types_wt),
         source = factor(source, levels = c("Genes", "Enhancers")))

write.csv(enrichR_df_wt, "../results/RNA_ATAC_integration/tables/KEGG_enrichment_wt.csv", row.names = FALSE)
# 
# enrichR_df_wt <- read.csv("../results/RNA_ATAC_integration/tables/KEGG_enrichment_wt.csv")
# cell_types_wt <- unique(enrichR_df_wt$cell_type)
# enrichR_df_wt <- enrichR_df_wt %>%
#   mutate(cell_type = factor(cell_type, levels = cell_types_wt),
#          source = factor(source, levels = c("Genes", "Enhancers")))


enrichR_df_ko <- bind_rows(enrichR_list_ko, .id = "cell_type") %>% 
  mutate(cell_type = factor(cell_type, levels = cell_types_ko),
         source = factor(source, levels = c("Genes", "Enhancers")))

write.csv(enrichR_df_ko, "../results/RNA_ATAC_integration/tables/KEGG_enrichment_ko.csv", row.names = FALSE)

# enrichR_df_ko <- read.csv("../results/RNA_ATAC_integration/tables/KEGG_enrichment_ko.csv") 
# cell_types_ko <- unique(enrichR_df_ko$cell_type)
# enrichR_df_ko <- enrichR_df_ko%>%
#   mutate(cell_type = factor(cell_type, levels = cell_types_ko),
#          source = factor(source, levels = c("Genes", "Enhancers")))


top_terms_wt <- enrichR_df_wt %>%
  filter(status == "Activated") %>% 
  group_by(cell_type, source) %>%
  slice_max(order_by = log10_Pval, n = 10) %>% 
  pull(Term) %>% 
  unique()

top_terms_ko <- enrichR_df_ko %>%
  filter(status == "Activated") %>% 
  group_by(cell_type, source) %>%
  slice_max(order_by = log10_Pval, n = 10) %>%
  pull(Term) %>% 
  unique()

new_names_wt <- c("Artery", "Capillary artery", "Capillary", "Capillary vein", "Large vein")
names(new_names_wt) <- levels(enrichR_df_wt$cell_type)

new_names_ko <- c("Artery KO", "Capillary artery KO", "Capillary KO", "Capillary vein KO", "Large vein KO")
names(new_names_ko) <- levels(enrichR_df_ko$cell_type)


data_plot_wt <- enrichR_df_wt %>% 
  mutate(cell_type = recode(as.character(cell_type), !!!new_names_wt)) %>% 
  filter(Term %in% top_terms_wt,
         status == "Activated") %>% 
  mutate(cell_type = factor(cell_type, levels = new_names_wt)) %>% 
  arrange(source, cell_type, desc(log10_Pval)) %>% 
  mutate(Term = fct_rev(fct_inorder(Term))) 

data_plot_ko <- enrichR_df_ko %>% 
  mutate(cell_type = recode(as.character(cell_type), !!!new_names_ko)) %>% 
  filter(Term %in% top_terms_ko,
         status == "Activated") %>% 
  mutate(cell_type = factor(cell_type, levels = new_names_ko)) %>% 
  arrange(source, cell_type, desc(log10_Pval)) %>% 
  mutate(Term = fct_rev(fct_inorder(Term))) 


EC_colors <- c("#247567", 
               "#E2705B", 
               "#6C7A99", 
               "#B4698B", 
               "#D8A21F")

plot_dot_kegg <- function(data) {
  ggplot(data, aes(x = cell_type, y = Term, color = log10_Pval, size = Count)) +
    geom_point() +
    scale_color_viridis_c(breaks = seq(floor(min(data$log10_Pval)), ceiling(max(data$log10_Pval)), length.out = 5),
                          labels = function(x) round(x, 0),name = "-log10\nP-value") +
    theme_bw() +
    ggtitle("KEGG enrichment") +
    theme(axis.text.y = element_text(size = 9, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "bold", colour = EC_colors),
          legend.key.height = unit(12, "pt"),
          legend.title = element_text(size = 10, face = "bold", margin = margin(t = 2, b = 4)),
          legend.text = element_text(size = 9),
          axis.ticks.x = element_line(linewidth = 0),
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 9.5, face = "bold", hjust = 0.5, vjust = 0.5, colour = "black")) +
    scale_size_continuous(range = c(1, 3)) +
    facet_wrap(~source, scales = "free_x")
}


plot_dot_kegg(data_plot_wt)
plot_dot_kegg(data_plot_ko)

ggsave("../results/RNA_ATAC_integration/graphs/kegg_dotplot_wt.png", plot = plot_dot_kegg(data_plot_wt),
       width = 6, height = 7, dpi = 300)

ggsave("../results/RNA_ATAC_integration/graphs/kegg_dotplot_ko.png", plot = plot_dot_kegg(data_plot_ko),
       width = 6, height = 6, dpi = 300)


selected_terms <- read.csv("utils/selected_terms.csv") %>% pull("selected_pathways")


selected_to_plot <- enrichR_df_ko %>% 
  mutate(cell_type = recode(as.character(cell_type), !!!new_names_ko)) %>% 
  filter(Term %in% selected_terms,
         status == "Activated") %>% 
  mutate(cell_type = factor(cell_type, levels = new_names_ko),
         Term = fct_rev(fct_inorder(Term)))

plot_dot_kegg(selected_to_plot)


ggsave("../results/RNA_ATAC_integration/graphs/kegg_dotplot_ko_selected.png", plot = plot_dot_kegg(selected_to_plot),
       width = 5.2, height = 4.7, dpi = 300)


unique_terms_wt <- enrichR_df_wt %>%
  filter(status == "Activated",
         source == "Genes") %>% 
  select(cell_type, Term)  %>%
  group_by(Term) %>% 
  summarise(n_celltypes = n()) %>% 
  filter(n_celltypes == 1) %>% pull(Term)


specific_terms_plot <-  enrichR_df_wt %>% 
  mutate(cell_type = recode(as.character(cell_type), !!!new_names_wt)) %>% 
  filter(Term %in% unique_terms_wt,
         status == "Activated") %>% 
  mutate(cell_type = factor(cell_type, levels = new_names_wt)) %>% 
  arrange(source, cell_type, desc(log10_Pval)) %>% 
  mutate(Term = fct_rev(fct_inorder(Term)))

plot_dot_kegg(specific_terms_plot)


ggsave("../results/RNA_ATAC_integration/graphs/kegg_dotplot_wt_specific.png", plot = plot_dot_kegg(specific_terms_plot),
       width = 6.2, height = 7.5, dpi = 300)


# === Log R session info for reproducibility ===
sink("../results/logs/sessioninfo_16_KEGG_enrichment.txt")
sessionInfo()
sink()




