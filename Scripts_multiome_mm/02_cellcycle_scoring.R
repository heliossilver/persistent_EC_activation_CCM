# Libraries
pacman::p_load("Seurat", "dplyr", "ggplot2", "cowplot", "glue", "purrr")

# Create output folders if they don't exist
dir.create("../data/databases_human_mouse_integration", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/seurat_objects/filtered_ccscore", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/plots/ccscore", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/plots/graphs_endo", recursive = TRUE, showWarnings = FALSE)

# Load or download JAX ortholog table
jax_file <- "../data/databases_human_mouse_integration/jax_ms_hu_orthologs.csv"
if (!file.exists(jax_file)) {
  mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep = "\t")
  write.csv(mouse_human_genes, jax_file, row.names = FALSE)
} else {
  mouse_human_genes <- read.csv(jax_file)
}

# Cell cycle gene conversion function
convert_human_to_mouse <- function(gene_list){
  human_filtered <- mouse_human_genes %>%
    filter(Symbol %in% gene_list, Common.Organism.Name == "human")
  mouse_filtered <- mouse_human_genes %>%
    filter(Common.Organism.Name == "mouse, laboratory")
  merged <- left_join(human_filtered, mouse_filtered, by = "DB.Class.Key") %>%
    select(Symbol.x, Symbol.y) %>%
    rename(gene_hu = Symbol.x, gene_ms = Symbol.y)
  unique(merged$gene_ms)
}

# Get Seurat cell cycle gene sets and convert
s.genes.h <- cc.genes.updated.2019$s.genes
g2m.genes.h <- cc.genes.updated.2019$g2m.genes

s.genes.m <- convert_human_to_mouse(s.genes.h)
g2m.genes.m <- convert_human_to_mouse(g2m.genes.h)

# Load filtered Seurat objects (with consistent capitalization)
sample_names <- c("Flox1", "Flox2", "Flox3", "KO1", "KO2", "KO3")
endo.list <- setNames(lapply(sample_names, function(name) {
  readRDS(glue("../results/seurat_objects/filtered/mBMEC_{name}_multiome_filtered.rds"))
}), sample_names)

# Score cell cycle
endo.list <- lapply(endo.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x))
  x <- RunPCA(x)
  x <- FindNeighbors(x, dims = 1:10)
  x <- FindClusters(x, resolution = 0.5)
  x <- RunUMAP(x, dims = 1:10, verbose = FALSE)
  x <- CellCycleScoring(x, s.features = s.genes.m, g2m.features = g2m.genes.m)
  return(x)
})

# Save scored objects with consistent naming
walk2(names(endo.list), endo.list, ~ {
  saveRDS(.y, glue("../results/seurat_objects/filtered_ccscore/{.x}_filtered_ccscore.rds"))
})

# Reload scored objects (if needed)
endo.list <- setNames(lapply(sample_names, function(name) {
  readRDS(glue("../results/seurat_objects/filtered_ccscore/{name}_filtered_ccscore.rds"))
}), sample_names)

# Merge objects with cell IDs
endo <- merge(x = endo.list[[1]],
              y = endo.list[-1],
              add.cell.ids = sample_names)

 
# Run PCA on cell cycle genes
endo <- ScaleData(endo, features = rownames(endo)) %>% #c(s.genes.m, g2m.genes.m)) %>% 
  RunPCA(features = c(s.genes.m, g2m.genes.m), approx = FALSE)

# Generate and save PCA plot split by phase
cell_cycle_plot <- DimPlot(object = endo,
                           reduction = "pca",
                           group.by = "Phase") 

cell_cycle_plot_split_by_phase <- DimPlot(object = endo,
                                          reduction = "pca",
                                          group.by = "Phase",
                                          split.by = "orig.ident")

ggsave(filename = "../results/plots/graphs_endo/cell_cycle.png",
       plot = cell_cycle_plot,
       device = "png",
       width = 4,
       height = 4,
       units = "in",
       dpi = 300)

ggsave(filename = "../results/plots/graphs_endo/cell_cycle_split_by_phase.png",
       plot = cell_cycle_plot_split_by_phase,
       device = "png",
       width = 6,
       height = 4,
       units = "in",
       dpi = 300)

# ==== Reproducibility log ====
sink("../logs/sessioninfo_02_cellcycle_scoring.txt")
sessionInfo()
sink()
