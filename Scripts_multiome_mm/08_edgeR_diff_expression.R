# Script 07: Differential expression with edgeR from pseudobulk RNA
pacman::p_load("edgeR", "limma", "ggplot2", "dplyr", "glue", "readr", 
               "ggplot2",  "stringr", "tidyr", "tibble", "purrr", 
               "ggtext", "EnhancedVolcano", "ComplexUpset")

# Output folders----
deg_folder <- "../results/RNA/edgeR/deg_tables"
deg_difexp_folder <- "../results/RNA/edgeR/deg_difexp_tables"
rds_folder <- "../results/RNA/rds"

dir.create(deg_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_difexp_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_folder, showWarnings = FALSE)

# Load counts----
cts <- read.csv("../results/RNA/data/counts_rna_all.csv", row.names = 1)
colnames(cts) <- gsub("\\.", "_", colnames(cts))
colnames(cts) <- paste0(colnames(cts), "_RNA")

# Calculate CPM and log2CPM
cpm <- edgeR::cpm(cts)
log2cpm <- log2(cpm + 1)
avg_logCPM <- rowMeans(log2cpm)

# hist(avg_logCPM, breaks=50, main="Histogram of Average logCPM",
#      xlab="Average logCPM", col="lightblue", border="black")
# abline(v=1, col="red", lwd=2, lty=2)
# 
# 

# Save CPM tables
write_csv(rownames_to_column(as.data.frame(cpm), var = "gene"), "../results/RNA/data/cpm.csv")
write_csv(rownames_to_column(as.data.frame(log2cpm), var = "gene"), "../results/RNA/data/log2cpm.csv")

# Design
ctsGroups <- c("AF", "AF", "AF", "AK", "AK", "AK", 
               "cAF", "cAF", "cAF", "cAK", "cAK", "cAK",
               "cF", "cF", "cF", "cK", "cK", "cK",
               "cVF", "cVF", "cVF", "cVK", "cVK", "cVK",
               "LVF", "LVF", "LVF", "LVK", "LVK", "LVK")

# edgeR pipeline
d <- DGEList(counts = cts, group = factor(ctsGroups))
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)

# Filter genes
keep <- avg_logCPM > 0.5
d <- d[keep, , keep.lib.sizes = FALSE]

# Save filtered CPM tables
write_csv(filter(rownames_to_column(as.data.frame(cpm), var = "gene"), gene %in% rownames(d)), "../results/RNA/data/cpm_filtered.csv")
write_csv(filter(rownames_to_column(as.data.frame(log2cpm), var = "gene"), gene %in% rownames(d)), "../results/RNA/data/log2cpm_filtered.csv")

#mds
# par(mar = c(4, 4, 3, 4), xpd = TRUE)  # Add space on the right with mar[4]
# plotMDS(d, col = as.numeric(d$samples$group), pch = 20)
# legend("topright", inset = c(-0.27, 0), 
#        legend = as.character(unique(d$samples$group)), 
#        col = 1:10, pch = 20, xpd = TRUE)
# Design matrix and model
design <- model.matrix(~ 0 + d$samples$group)
colnames(design) <- levels(d$samples$group)
d <- estimateDisp(d, design)
fit <- glmFit(d, design)
# plotBCV(d)

# Contrasts
contrast_mat <- makeContrasts(Artery_vs_rest        = AF  - (cAF + cVF + cF + LVF)/4,
                              CapArt_vs_rest        = cAF - (AF + cVF + cF + LVF)/4,
                              Cap_vs_rest           = cF  - (AF + cAF + cVF + LVF)/4,
                              CapVein_vs_rest       = cVF - (AF + cAF + cF + LVF)/4,
                              LargeVein_vs_rest     = LVF - (AF + cAF + cVF + cF)/4,
                              Artery_KO_vs_Flox     = AK  - AF,
                              CapArt_KO_vs_Flox     = cAK - cAF,
                              Cap_KO_vs_Flox        = cK  - cF,
                              CapVein_KO_vs_Flox    = cVK - cVF,
                              LargeVein_KO_vs_Flox  = LVK - LVF,
                              levels = design)

# Run tests
dge.l <- list()
contrasts.l <- list()

for (contrast_name in colnames(contrast_mat)) {
  lrt <- glmLRT(fit, contrast = contrast_mat[, contrast_name])
  result <- topTags(lrt, n = Inf, adjust.method = "fdr")$table %>%
    rownames_to_column("gene") %>%
    rename_with(~ paste0(., "_", contrast_name), -gene)
  
  dge.l[[contrast_name]] <- result
  contrasts.l[[contrast_name]] <- lrt
  
  write_csv(result, file = file.path(deg_folder, paste0(contrast_name, ".csv")))
}

# Save full lists
write_csv(reduce(dge.l, full_join, by = "gene"), file.path(deg_folder, "allcelltypes.csv"))
saveRDS(dge.l, file = file.path(rds_folder, "dge_celltype_list.rds"))
saveRDS(contrasts.l, file = file.path(rds_folder, "contrasts_list.rds"))


# Annotate up/down
for (name in names(dge.l)) {
  df <- dge.l[[name]]
  logFC_col <- paste0("logFC_", name)
  FDR_col <- paste0("FDR_", name)

  df <- df %>%
    mutate(diffexpressed = case_when(
      !!sym(FDR_col) <= 0.05 & !!sym(logFC_col) >= 0.25 ~ "Up",
      !!sym(FDR_col) <= 0.05 & !!sym(logFC_col) <= -0.25 ~ "Down",
      TRUE ~ "No"
    )) %>%
    arrange(desc(!!sym(logFC_col)))

  write_csv(df, file = file.path(deg_difexp_folder, paste0(name, "_diffexp_full.csv")))
  write_csv(filter(df, diffexpressed != "No"), file = file.path(deg_difexp_folder, paste0(name, "_diffexp.csv")))
}

# =========================
# Volcano plots (all contrasts)
# =========================


volcano_dir <- "../results/plots/graphs_RNA"
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

# Display name mapping
full_names <- c("Artery" = "Artery",
  "CapArt" = "Capillary artery",
  "Cap" = "Capillary",
  "CapVein" = "Capillary vein",
  "LargeVein" = "Large vein")

# Loop and plot
for (name in names(dge.l)) {
  df <- dge.l[[name]]

  logFC_col <- paste0("logFC_", name)
  FDR_col <- paste0("FDR_", name)

  df <- df %>%
    column_to_rownames("gene") %>%
    arrange(desc(!!sym(logFC_col))) %>%
    mutate(diffexpressed = case_when(!!sym(FDR_col) <= 0.05 & !!sym(logFC_col) >= 0.25 ~ "Up",
                                     !!sym(FDR_col) <= 0.05 & !!sym(logFC_col) <= -0.25 ~ "Down",
                                     TRUE ~ "No"))

  summary_data <- df %>%
    group_by(diffexpressed) %>%
    summarise(n = n(), .groups = "drop")

  abbrev <- stringr::str_extract(name, "^[A-Za-z]+")
  full_name <- full_names[abbrev]
  subtitle <- glue::glue("{summary_data$n[summary_data$diffexpressed == 'Down']} downregulated and {summary_data$n[summary_data$diffexpressed == 'Up']} upregulated")
  comparison_type <- if (stringr::str_detect(name, "KO_vs_Flox")) "DEG KO vs Flox" else "DEGs"
  
  volcano_p <- EnhancedVolcano(df,
                               lab = rownames(df),
                               x = logFC_col,
                               y = FDR_col,
                               title = glue::glue("{full_name} {comparison_type}"),
                               subtitle = subtitle,
                               subtitleLabSize = 15,
                               pCutoff = 0.05,
                               FCcutoff = 0.25,
                               legendLabSize = 12,
                               axisLabSize = 13)

  ggsave(filename = file.path(volcano_dir, paste0(name, "_volcano.png")),
         plot = volcano_p, width = 6, height = 5, dpi = 300)
}

# =========================
# DEG summary barplots
# =========================


summary_dir <- "../results/plots/graphs_endo"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# Load differential expression files
deg_files <- list.files("../results/RNA/edgeR/deg_difexp_tables/", pattern = "_diffexp.csv", full.names = TRUE)
samples <- basename(deg_files) %>% str_remove("_diffexp.csv")

# Build summary table
gene_df_l <- map2(deg_files, samples, ~ {
  read.csv(.x) %>%
    rename(logFC = matches("^log[_]?FC")) %>%
    filter(abs(logFC) > 0.5) %>% 
    select(gene, diffexpressed) %>%
    mutate(cell_type = str_replace(.y, "_vs.*", ""))
})

gene_df <- bind_rows(gene_df_l)

gene_summ <- gene_df %>%
  count(cell_type, diffexpressed, name = "num_genes") %>%
  mutate(condition = if_else(str_detect(cell_type, "KO"), "KO", "Flox"), 
         cell_type = str_replace_all(cell_type, 
                                     c("_KO" = "",
                                       "CapArt" = "Capillary artery", 
                                       "CapVein|CapVen" = "Capillary vein",
                                       "LargeVein" = "Large vein", 
                                       "\\bCap\\b" = "Capillary")),
         diffexpressed = recode(diffexpressed,
                                "Up" = "Up\nregulated",
                                "Down" = "Down\nregulated"))

write.csv(gene_summ, "../results/RNA/data/summary_DEGs.csv", row.names = FALSE)

# Color palette
EC_colors <- rev(c("#247567", "#E2705B", "#6C7A99", "#B4698B", "#D8A21F"))

# -------- FLOX
gene_summ_Fl <- filter(gene_summ, condition == "Flox")
gene_summ_Fl$cell_type <- factor(gene_summ_Fl$cell_type, levels = rev(c("Artery", "Capillary artery",
                                                                        "Capillary", "Capillary vein",
                                                                        "Large vein")))
gene_summ_Fl$diffexpressed <- factor(gene_summ_Fl$diffexpressed, levels = rev(unique(gene_summ_Fl$diffexpressed)))

gene_fl_p <- ggplot(gene_summ_Fl, aes(y = cell_type)) +
  geom_segment(aes(x = 0,
                   xend = ifelse(diffexpressed == "Down\nregulated", -num_genes, num_genes),
                   yend = cell_type,
                   color = diffexpressed), linewidth = 1) +
  geom_point(aes(x = ifelse(diffexpressed == "Down\nregulated", -num_genes, num_genes),
                 color = diffexpressed), size = 3) +
  geom_text(aes(x = ifelse(diffexpressed == "Down\nregulated", -num_genes, num_genes),
                label = abs(num_genes),
                hjust = ifelse(diffexpressed == "Down\nregulated", 1.25, -0.39)),
            size = 3.5) +
  scale_color_manual(values = c("Down\nregulated" = "dodgerblue", "Up\nregulated" = "#EE6F5D")) +
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_text(size = 11, margin = margin(t = 8, unit = "pt")),
         axis.text.y = element_text(size = 10, hjust = 1, face = "bold", colour = EC_colors),
         legend.text = element_text(size = 9.5, face = "bold", margin = margin(b = 5)),
         axis.title.y = element_blank(),
         plot.title = ggtext::element_markdown(hjust = 0.5),
         plot.subtitle = element_text(hjust = 0.5, size = 9)) +
  labs(x = "Number of significant DEGs", y = "", color = "",
       title = "<b>BEC subtype marker genes</b><br><b>(</b><i><b>Pdcd10</b></i><sup><b>fl/fl</b></sup><b>)</b>",
       subtitle = "Up/down-regulated DEGs |FC| > 0.5") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25) +
  scale_x_continuous(breaks = seq(-1250, 1250, by = 500), labels = function(x) abs(x), limits = c(-1500, 1700)) 


ggsave(file.path(summary_dir, "lollipop_genes_Fl.png"), gene_fl_p, width = 5, height = 3, dpi = 500)

# -------- KO
gene_summ_KO <- filter(gene_summ, condition == "KO")
gene_summ_KO$cell_type <- paste(gene_summ_KO$cell_type, "KO")
gene_summ_KO$cell_type <- factor(gene_summ_KO$cell_type, levels = rev(paste(c("Artery", "Capillary artery",
                                                                                           "Capillary", "Capillary vein",
                                                                                           "Large vein"), "KO")))
gene_summ_KO$diffexpressed <- factor(gene_summ_KO$diffexpressed, levels = rev(unique(gene_summ_KO$diffexpressed)))

gene_ko_p <- ggplot(gene_summ_KO, aes(y = cell_type)) +
  geom_segment(aes(x = 0,
                   xend = ifelse(diffexpressed == "Down\nregulated", -num_genes, num_genes),
                   yend = cell_type,
                   color = diffexpressed), linewidth = 1) +
  geom_point(aes(x = ifelse(diffexpressed == "Down\nregulated", -num_genes, num_genes),
                 color = diffexpressed), size = 3) +
  geom_text(aes(x = ifelse(diffexpressed == "Down\nregulated", -num_genes, num_genes),
                label = abs(num_genes),
                hjust = ifelse(diffexpressed == "Down\nregulated", 1.25, -0.39)),
            size = 3.5) +
  scale_color_manual(values = c("Down\nregulated" = "dodgerblue", "Up\nregulated" = "#EE6F5D")) +
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_text(size = 11, margin = margin(t = 8, unit = "pt")),
        axis.text.y = element_text(size = 10, hjust = 1, face = "bold", colour = EC_colors),
        legend.text = element_text(size = 9.5, face = "bold", margin = margin(b = 5)),
        axis.title.y = element_blank(),
        plot.title = ggtext::element_markdown(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9)) +
  labs(x = "Number of significant DEGs", y = "", color = "",
       title = "<b>DEG per BEC subtype</b><br><b>(</b><i><b>Pdcd10</b></i><sup><b>BECKO</b></sup> vs <i><b>Pdcd10</b></i><sup><b>fl/fl</b></sup><b>)</b>",
       subtitle = "Up/down-regulated DEGs |FC| > 0.5") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25) +
  scale_x_continuous(breaks =  seq(-3000, 3000, by = 1500), labels = function(x) abs(x), limits = c(-3000, 4000)) 


ggsave(file.path(summary_dir, "lollipop_genes_KO.png"), gene_ko_p, width = 4.5, height = 3, dpi = 500)

# =========================
# Upset plots KO
# =========================

# Load differential expression files
deg_files <- list.files("../results/RNA/edgeR/deg_difexp_tables/", pattern = "_vs_Flox_diffexp.csv", full.names = TRUE)
samples <- basename(deg_files) %>% str_remove("_vs_Flox_diffexp.csv")

# Build summary table
gene_df_l <- map2(deg_files, samples, ~ {
  read.csv(.x) %>%
    filter(diffexpressed == "Up") %>% pull("gene") %>% unique()
})

gene_df_l <- gene_df_l[c(1, 3, 2, 4, 5)]

names(gene_df_l) <- c("Artery", "Capillary artery","Capillary","Capillary vein","Large vein")

from_list <- function(list_data) {
  members = unique(unlist(list_data))
  data.frame(
    lapply(list_data, function(set) members %in% set),
    row.names=members,
    check.names=FALSE
  )
}


set_order <- rev(c("Artery", "Capillary artery","Capillary","Capillary vein","Large vein"))

gene_up_mtx <- from_list(gene_df_l)

upset_gene_up <- upset(gene_up_mtx, set_order, n_intersections = 10, 
                       keep_empty_groups = T, sort_sets = F, name = NULL,
                       base_annotations = list('Intersection size'= (intersection_size(text = list(size = 4.5, vjust = -0.2),
                                                                                       bar_number_threshold = 1,
                                                                                       fill='grey80') ) +
                                                 labs(title = "<b>Intersection of upregulated DEGs</b><br>
                                                <b>(</b><i><b>PDCD10</b></i><sup><b>BECKO</b></sup> 
                                                     <b> vs </b> <i><b>PDCD10</b></i><sup><b>fl/fl</b></sup><b>)</b>") + 
                                                 scale_y_continuous(limits = c(0, 600)) +
                                                 theme(text = element_text(size = 12), 
                                                       axis.text.y = element_text(size = 14),
                                                       axis.title.y = element_text(size = 15, vjust = -20),
                                                       plot.title = ggtext::element_markdown(hjust = 0.5, size = 12))),
                       height_ratio = 0.3, set_sizes = F) + theme(axis.text.y = element_text(size = 12, face = "bold"))


ggsave(filename = "../results/plots/graphs_endo/upset_gene_up_KO.png", plot = upset_gene_up, width = 7, height = 5, bg = "white")


# Optionally save session info
sink("../logs/sessioninfo_08_edgeR_diff_expression.txt")
sessionInfo()
sink()