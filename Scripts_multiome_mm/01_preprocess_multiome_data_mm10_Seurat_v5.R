#Libraries

pacman::p_load("Seurat", "Signac", "harmony", 
               "dplyr", "ggplot2", "rtracklayer",
               "DoubletFinder", "EnsDb.Mmusculus.v79")

options(future.globals.maxSize = 2000 * 1024^2)  # 2000 MiB (2 GiB)

# Create output folders if they don't exist
dir.create("../results/seurat_objects/unfiltered", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/seurat_objects/filtered", recursive = TRUE, showWarnings = FALSE)

# Load gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"


list.files("../multiome_h5_fragments/")

# List of sample names
samples <- list("Flox1", "Flox2", "Flox3", "KO1", "KO2", "KO3")

# Process each sample
for (loop_var in samples) {
  tryCatch({
    message(glue::glue("⏳ Processing {loop_var}..."))
    
    # 1. Load raw counts (RNA and ATAC)
    data <- Read10X_h5(sprintf('../multiome_h5_fragments/%s_filtered_feature_bc_matrix.h5',loop_var))
    rna_counts <- data$`Gene Expression`
    atac_counts <- data$Peaks
    
    # 2. Create RNA Seurat object
    obj <- CreateSeuratObject(counts = rna_counts, project = sprintf('%s',loop_var))
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
    
    # 3. Filter ATAC peaks (standard chromosomes only)
    grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]
    
    # 4. Add ATAC assay
    frag.file <- sprintf('../multiome_h5_fragments/%s_atac_fragments.tsv.gz',loop_var)
    chrom_assay <- CreateChromatinAssay(counts = atac_counts,
                                        sep = c(":", "-"), 
                                        genome = 'mm10',
                                        fragments = frag.file, 
                                        min.cells = 10, 
                                        annotation = annotations)
    
    obj[["ATAC"]] <- chrom_assay
    obj$orig.ident <- sprintf('%s',loop_var)
    
    # 5. Calculate QC metrics (TSS enrichment) and filter cells
    DefaultAssay(obj) <- "ATAC"
    obj <- TSSEnrichment(obj)
    DefaultAssay(obj) <- "RNA"
    
    # Save unfiltered object
    saveRDS(obj, file = sprintf('../results/seurat_objects/unfiltered/mBMEC_%s_multiome_unfiltered.rds', loop_var))
    
    # Apply basic RNA/ATAC QC filtering
    obj <- subset(x = obj, 
                  subset = nCount_ATAC > 500 & nCount_RNA < 50000 & nFeature_RNA > 200 & percent.mt < 20)
    
    # 6. Doublet removal
    nExp <- round(ncol(obj) * 0.1)
    
    set.seed("1234")
    obj <- doubletFinder(obj, pN = 0.25, pK = 0.01, nExp = nExp, PCs = 1:20, sct = TRUE)
    DF.name <- colnames(obj@meta.data)[grepl("DF.classification", colnames(obj@meta.data))]
    obj <- obj[, obj@meta.data[, DF.name] == "Singlet"]

    # Save filtered object
    saveRDS(obj, file = sprintf('../results/seurat_objects/filtered/mBMEC_%s_multiome_filtered.rds', loop_var))
    
    message(glue::glue("✅ Finished {loop_var}"))
  }, error = function(e) {
    message(glue::glue("❌ Error processing {loop_var}: {e$message}"))
  })
}

# ==== Reproducibility log ====
dir.create("../results/logs", showWarnings = FALSE)
sink("../logs/sessioninfo_01_preprocess_multiome_data_mm10_Seurat_v5.txt")
sessionInfo()
sink()


