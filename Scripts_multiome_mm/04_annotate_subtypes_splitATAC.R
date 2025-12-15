# Libraries
pacman::p_load("Seurat", "Signac", "vroom",  "tidyr", "scCustomize",
               "dplyr", "ggplot2", "RColorBrewer",  "tibble", 
               "dittoSeq", "glue", "ggtext", "clustree")

# Directories
dir.create("../results/plots/graphs_endo", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/ATAC/bedfiles_signac", recursive = TRUE, showWarnings = FALSE)

# Load filtered mapped object
endo_filt <- readRDS("../results/seurat_objects/merged/endo_filtered_mapped.rds")


# Re-cluster at high resolution for subtype dissection
endo_filt <- FindNeighbors(endo_filt, reduction = "integrated.cca", dims = 1:15, graph.name = "res") %>% 
  FindClusters(resolution = 5)


# UMAP
endo5_UMAP.p <- DimPlot(endo_filt, reduction = "umap", pt.size = 1.2, alpha = 0.7, label.size = 5,
                        label = TRUE, repel = TRUE, label.box = TRUE, group.by = "RNA_snn_res.5") +
  labs(title = "endo_Res_5") +
  theme(plot.title = element_text(hjust = 0.5, size = 16), legend.position = "none")

# Cluster composition plot
bar_data <- dittoBarPlot(endo_filt, var = "predicted.EC_type", group.by = "RNA_snn_res.5", sub = "5", 
                         retain.factor.levels = TRUE, do.hover = TRUE, data.out = TRUE)
endo_filt3_table <- bar_data$data

# Annotate based on majority label per cluster
top_percentages <- endo_filt3_table %>%
  group_by(grouping) %>%
  slice_max(order_by = percent, n = 1) %>%
  mutate(final_label = case_when(percent > 0.50 ~ as.character(label), TRUE ~ "unknown")) 


endo_filt_subcluster_table <- bar_data_subcluster$data

top_percentages_final <- top_percentages %>%
  mutate(final_label = case_when(grouping == "14" ~ "capillary vein",
                                 grouping == "36" ~ "capillary vein",
                                 grouping == "39" ~ "artery",
                                 grouping == "50" ~ "artery",
                                 grouping == "53" ~ "unknown",
                                 grouping == "54" ~ "unknown",
                                 grouping == "55" ~ "unknown",
                                 TRUE ~ final_label)) %>% 
  mutate(final_label = stringr::str_replace_all(final_label, c("arterial" = "artery", "venous" = "vein")))

# Rename clusters
new_clusters <- setNames(top_percentages_final$final_label, top_percentages_final$grouping)
endo_filt <- RenameIdents(endo_filt, new_clusters)
endo_filt[["cell_type"]] <- Idents(endo_filt)

# UMAPs: annotated + splits
endoannotated_UMAP.p <- DimPlot(endo_filt, reduction = "umap", pt.size = 0.8, group.by = "cell_type") +
  labs(title = "Brain endothelial cells") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), legend.text = element_text(size = 18))

endoannotated_UMAP_split.p <- DimPlot(endo_filt, reduction = "umap", group.by = "cell_type", split.by = "condition", pt.size = 0.8, label.size = 5, repel = TRUE) +
  labs(title = "Brain endothelial cells") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), strip.text = element_text(size = 18), legend.text = element_text(size = 18))

# Remove CP, reprocess
endo <- subset(endo_filt, subset = cell_type %in% c("unknown", "choroid plexus"), invert = TRUE)
endo <- NormalizeData(endo) %>% 
  FindVariableFeatures() %>%
  ScaleData(features = rownames(endo)) %>%
  RunPCA()
ElbowPlot(endo)

endo <- endo %>%
  FindNeighbors(reduction = "integrated.cca", dims = 1:15) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:13, reduction = "integrated.cca")


# Set cell type order and palette
EC_colors <- c("#247567", "#E2705B", "#6C7A99", "#B4698B", "#D8A21F")
endo$cell_type <- factor(endo$cell_type, levels = c("artery", "capillary artery", "capillary", "capillary vein", "large vein"))
Idents(endo) <- "cell_type"

# Annotated UMAPs
endoannotated_UMAP.p <- DimPlot(endo, reduction = "umap", label = TRUE, label.box = TRUE, label.color = "white", repel = TRUE,
                                cols = EC_colors, group.by = "cell_type", label.size = 4.5) +
  labs(title = "Annotated Endothelial Cells") +
  theme(plot.title = element_text(hjust = 0.5, size = 17))
endoannotated_UMAP.p

ggsave("../results/plots/graphs_endo/endoannotated_UMAP_noCP.png", endoannotated_UMAP.p, width = 6, height = 5, dpi = 300)

# By condition
condition_umap <- DimPlot(endo, reduction = "umap", group.by = "condition") +
  labs(caption = paste0(ncol(endo), " nuclei")) +
  theme(text = element_text(size = 20))
ggsave("../results/plots/graphs_endo/endoannotated_condition_UMAP_noCP.png", condition_umap, width = 5, height = 5, dpi = 300)


## Subset by condition
flox <- subset(endo, subset = condition == "Flox")
ko <-  subset(endo, subset = condition == "KO")

# Dot plot of markers
DefaultAssay(endo) <- "RNA"
dotplot_end <- DotPlot(endo, features = c('Fbln5', 'Gkn3', 'Hey1', "Igfbp4","Bmx",
                                          'Tgfb2', "Stmn2", "Alpl",
                                          'Mfsd2a', 'Rgcc', 'Tfrc', "Car4", "Slc22a8", "Slc16a1",
                                          'Slc38a5', 'Vwf', 'Vcam1', "Ttyh2",  "Pam", "Irf1"),
                       cols = c("salmon1", "salmon1"), dot.scale = 4, split.by = "condition") +
  geom_vline(xintercept = c(5.5, 8.5, 10.5, 14.5, 20.5), linetype = "dashed", color = EC_colors, linewidth = 0.75) +
  coord_flip() +
  labs(title = "Endothelial cell markers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(),
        axis.title.y = element_blank(), title = element_text(size = 10), legend.title = element_text(size = 8))

ggsave("../results/plots/graphs_endo/dotplot_markersEC_noCP.png", dotplot_end, width = 4, height = 5, dpi = 300, bg = "white")

# Dot plot of markers Flox
DefaultAssay(flox) <- "RNA"
dotplot_flox <- DotPlot(flox, features = c('Fbln5', 'Gkn3', 'Hey1', "Igfbp4","Bmx",
                                          'Tgfb2', "Stmn2", "Alpl",
                                          'Mfsd2a', 'Rgcc', 'Tfrc', "Car4", "Slc22a8", "Slc16a1",
                                          'Slc38a5', 'Vwf', 'Vcam1', "Ttyh2",  "Pam", "Irf1"),
                       cols = c("lightblue1", "salmon1"), dot.scale = 4.5,
                       col.min = -1.8, col.max = 1.8) +
  geom_vline(xintercept = c(5.5, 8.5, 10.5, 14.5, 20.5), linetype = "dashed", color = EC_colors, linewidth = 0.75) +
  coord_flip() +
  labs(title = "<b>Endothelial cell markers</b><br><b>(</b><b><i>Pdcd10</b></i><sup><b>fl/fl</b></sup><b>)</b>") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        title = element_markdown(size = 12, face = "bold"), legend.title = element_text(size = 10), legend.text = element_text(size = 10))

ggsave("../results/plots/graphs_endo/dotplot_markersEC_Flox2.png", dotplot_flox, width = 4, height = 5, dpi = 300, bg = "white")

# Annotated UMAPs
flox_umap <- DimPlot(flox, reduction = "umap", label = TRUE, label.box = TRUE, label.color = "white",
                                cols = EC_colors, group.by = "cell_type", label.size = 5.5, repel = TRUE) +
  labs(title = "Annotated brain endothelial cells") +
  theme(plot.title = element_text(hjust = 0.5, size = 17), 
        legend.title = element_text(size = 12), legend.text = element_text(size = 12)) + NoLegend()

ggsave("../results/plots/graphs_endo/endoannotated_UMAP_flox.png", flox_umap, width = 5, height = 5, dpi = 300)
ggsave("../results/plots/graphs_endo/endoannotated_UMAP_flox.svg", flox_umap, width = 5, height = 5)
# Dot plot of markers KO
DefaultAssay(ko) <- "RNA"
dotplot_ko <- DotPlot(ko, features = c('Fbln5', 'Gkn3', 'Hey1', "Igfbp4","Bmx",
                                          'Tgfb2', "Stmn2", "Alpl",
                                          'Mfsd2a', 'Rgcc', 'Tfrc', "Car4", "Slc22a8", "Slc16a1",
                                          'Slc38a5', 'Vwf', 'Vcam1', "Ttyh2",  "Pam", "Irf1"),
                       cols = c("lightblue1", "salmon1"), dot.scale = 4.5, 
                      col.min = -1.8, col.max = 1.8) +
  geom_vline(xintercept = c(5.5, 8.5, 10.5, 14.5, 20.5), linetype = "dashed", color = EC_colors, linewidth = 0.75) +
  coord_flip() +
  labs(title = "<b>Endothelial cell markers</b><br><b>(</b><b><i>Pdcd10</b></i><sup><b>BECKO</b></sup><b>)</b>") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        title = element_markdown(size = 12, face = "bold"), legend.title = element_text(size = 10), legend.text = element_text(size = 10))

ggsave("../results/plots/graphs_endo/dotplot_markersEC_KO.png", dotplot_ko, width = 4, height = 5, dpi = 300, bg = "white")


# Cell type composition barplot (all samples)
endo$celltype_condition <- paste(endo$cell_type, endo$condition, sep = "_")

cell_condition_counts <- endo@meta.data %>%
  count(cell_type, orig.ident, name = "Count") %>%
  mutate(cell_type = factor(cell_type, levels = levels(endo$cell_type)))


summary_cell_type <- cell_condition_counts %>% 
  separate(orig.ident, into = c("condition", "replicate"), sep = -1) %>% 
  group_by(cell_type, condition) %>% 
  summarise(number_cells = sum(Count))
# 
# cell_type        condition number_cells
# <fct>            <chr>            <int>
#   1 artery           Flox               368
# 2 artery           KO                 172
# 3 capillary artery Flox               796
# 4 capillary artery KO                 446
# 5 capillary        Flox              5401
# 6 capillary        KO                3715
# 7 capillary vein   Flox              6110
# 8 capillary vein   KO                5189
# 9 large vein       Flox               572
# 10 large vein       KO                 469

sample_labels <- cell_condition_counts %>%
  group_by(orig.ident) %>%
  summarise(total_cells = sum(Count)) %>%
  mutate(label = paste0(orig.ident, "\n(", total_cells, ")")) %>%
  dplyr::select(orig.ident, label) %>%
  deframe()

# --- Plot: All replicates
barplot_all <- ggplot(cell_condition_counts, aes(x = orig.ident, y = Count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Cell Type Counts per Sample", x = "Sample", y = "Percentage",
       fill = "Cell Type") +
  scale_fill_manual(values = EC_colors) +
  scale_x_discrete(labels = sample_labels) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold", hjust = 0.5),
        plot.title = element_text(hjust = 0.5))

ggsave("../results/plots/graphs_endo/barplot_celltype_composition_all.png",
       barplot_all, width = 5, height = 4, dpi = 300)

# --- Plot: Flox only
barplot_flox <- ggplot(cell_condition_counts %>% dplyr::filter(grepl("Flox", orig.ident)),
                       aes(x = orig.ident, y = Count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Cell Type Counts per Sample", x = "Sample", y = "Percentage",
       fill = "Cell Type") +
  scale_fill_manual(values = EC_colors) +
  scale_x_discrete(labels = sample_labels) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold", hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(face = "bold"))

ggsave("../results/plots/graphs_endo/barplot_celltype_composition_flox.png",
       barplot_flox, width = 4, height = 3, dpi = 300)


# Split ATAC fragments
SplitFragments(endo,
               assay = "ATAC",
               group.by = "celltype_condition",
               outdir = "../results/ATAC/bedfiles_signac")


# Save final annotated object without CP
saveRDS(endo, "../results/seurat_objects/merged/endo_annotated_final.rds")

# ==== Reproducibility log ====
sink("../logs/sessioninfo_04_annotate_subtypes_splitATAC.txt")
sessionInfo()
sink()
