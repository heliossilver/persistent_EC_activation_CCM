
## ----Libraries-------------------------------------------------------------------------------------------------------------------------
pacman::p_load("Seurat", "Signac", "harmony", "dplyr", "ggplot2", "rtracklayer",
               "DoubletFinder", "EnsDb.Mmusculus.v79", "dittoSeq", "glue", "cowplot")

## ----Directories------------------------------------------------------------------------------------------------------------------------
dir.create("../results/seurat_objects/merged", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/plots/graphs_endo", recursive = TRUE, showWarnings = FALSE)

## ----Load unfiltered Seurat objects-----------------------------------------------------------------------------------------------------
sample_names <- c("Flox1", "Flox2", "Flox3", "KO1", "KO2", "KO3")

unfiltered <- setNames(lapply(sample_names, function(name) {
  readRDS(glue("../results/seurat_objects/unfiltered/mBMEC_{name}_multiome_unfiltered.rds"))
}), sample_names)

# Add condition metadata
for (name in names(unfiltered)) {
  unfiltered[[name]]$condition <- ifelse(grepl("Flox", name), "Flox", "KO")
}

# Merge unfiltered
unfilt <- merge(x = unfiltered[[1]], y = unfiltered[-1], add.cell.ids = sample_names)

saveRDS(unfilt, "../results/seurat_objects/merged/ccm_6samples_unfiltered_merged.rds")

## ----QC Violin Plots---------------------------------------------------------------------------------------------------------------------
vln_RNA.p   <- VlnPlot(unfilt, group.by = 'orig.ident', features = "nCount_RNA", log = TRUE, pt.size = 0)
vln_feat.p  <- VlnPlot(unfilt, group.by = 'orig.ident', features = "nFeature_RNA", log = TRUE, pt.size = 0)
vln_count.p <- VlnPlot(unfilt, group.by = 'orig.ident', features = "nCount_ATAC", log = TRUE, pt.size = 0)
vln_TSS.p   <- VlnPlot(unfilt, group.by = 'orig.ident', features = "TSS.enrichment", y.max = 20, pt.size = 0)
vln_mt.p    <- VlnPlot(unfilt, group.by = 'orig.ident', features = "percent.mt", pt.size = 0)


combined_plot <- plot_grid(vln_RNA.p, vln_feat.p, vln_mt.p, vln_count.p, vln_TSS.p, ncol = 3)

ggsave("../results/plots/graphs_endo/qc_violin_combined.png", combined_plot, width = 12, height = 8, dpi = 300)

## ----Load filtered Seurat objects with cell cycle scores----------------------------------------------------------------------------------
filtered <- setNames(lapply(sample_names, function(name) {
  readRDS(glue("../results/seurat_objects/filtered_ccscore/{name}_filtered_ccscore.rds"))
}), sample_names)

# Add condition metadata
for (name in names(filtered)) {
  filtered[[name]]$condition <- ifelse(grepl("Flox", name), "Flox", "KO")
}
# Merge filtered
endo <- merge(x = filtered[[1]], y = filtered[-1], add.cell.ids = sample_names)

saveRDS(endo, "../results/seurat_objects/merged/ccm_6samples_filtered_merged_ccscore.rds")

## ----Integration--------------------------------------------------------------------------------------------------------------------------
endo <- NormalizeData(endo) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA()

endo <- IntegrateLayers(object = endo,
                        method = CCAIntegration,
                        orig.reduction = "pca",
                        new.reduction = "integrated.cca",
                        verbose = FALSE)

endo[["RNA"]] <- JoinLayers(endo[["RNA"]])
ElbowPlot(endo)

endo <- FindNeighbors(endo, reduction = "integrated.cca", dims = 1:15)
endo <- FindClusters(endo, resolution = 0.1)
endo <- RunUMAP(endo, dims = 1:30, reduction = "integrated.cca")

saveRDS(endo, "../results/seurat_objects/merged/endo_integrated.rds")
endo <- readRDS("../results/seurat_objects/merged/endo_integrated.rds")
## ----Filter Non-Endothelial Cells---------------------------------------------------------------------------------------------------------
broad_markers_plot <- FeaturePlot(object = endo, reduction = "umap",
  features = c("Cdh5", "Pecam1", 
               "Acta2", 
               "Gfap", 
               "Col1a1", 
               "Ptprc",
               "Hbb-bs", 
               "Pdgfrb", "Myh11", "Vtn"),
  # c("Cdh5", "Pecam1", "Acta2", "Gfap", "Col1a1", "Ptprc"),
  repel = TRUE, ncol = 4)

broad_markers_vln <- VlnPlot(endo, features = c("Cdh5", "Pecam1", "Cldn5",
                                                "Acta2",  "Myh11",
                                                "Ptprc",
                                                "Hbb-bs", 
                                                "Pdgfrb", "Atp13a5"), pt.size = 0.1, ncol = 3)


ggsave("../results/plots/graphs_endo/broad_markers.png", broad_markers_plot, width = 10, height = 6, dpi = 300)
ggsave("../results/plots/graphs_endo/broad_markers_vln.png", broad_markers_vln, width = 12, height = 7, dpi = 300)
# Remove non-EC cells

endo_filt <- subset(endo, idents = setdiff(levels(endo), c("3", "4", "5")))

## ----Re-cluster after filtering-----------------------------------------------------------------------------------------------------------
endo_filt <- NormalizeData(endo_filt) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

endo_filt <- FindNeighbors(endo_filt, reduction = "integrated.cca", dims = 1:15)
endo_filt <- FindClusters(endo_filt,  resolution = 0.2)
endo_filt <- RunUMAP(endo_filt, dims = 1:15, reduction = "integrated.cca")

## ----Map to EC Reference (Kalucka 2022)---------------------------------------------------------------------------------------------------
source("Reference_processing.R")  # loads reference object

reference <- readRDS("../results/seurat_objects/reference_mapping_obj.rds")

DimPlot(reference, group.by = "Cluster", label = TRUE, repel = TRUE, label.box = TRUE)

anchors <- FindTransferAnchors(reference = reference,
                               query = endo_filt,
                               normalization.method = "LogNormalize",
                               reference.reduction = "pca",
                               dims = 1:50)

endo_filt <- MapQuery(anchorset = anchors,
  query = endo_filt,
  reference = reference,
  refdata = list(EC_type = "Cluster"))

saveRDS(endo_filt, "../results/seurat_objects/merged/endo_filtered_mapped.rds")

## ----Plot Mapped EC Types---------------------------------------------------------------------------------------------------------------
endo_EC_typecond <- DimPlot(endo_filt, reduction = "umap",
  split.by = "condition",
  label = TRUE, repel = TRUE, group.by = "predicted.EC_type",
  label.size = 8)

ggsave("../results/plots/graphs_endo/EC_types_mapped_split.png", endo_EC_typecond, height = 12, width = 12, dpi = 500)

# ==== Reproducibility log ====
sink("../results/logs/sessioninfo_03_merge_integrate_map.txt")
sessionInfo()
sink()
