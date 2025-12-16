# Persistent-Activation-of-BECs-linked-to-nueroinflammation-and-thrombosis-in-CCMs
This repository contains the scripts used to process data and generate figures and tables for:

**Gallego-Gutierrez et al.**  
*Persistent Activation of Endothelial Cells Is Linked to Neuroinflammation and Thrombosis in Cerebral Cavernous Malformation Disease*

The code is organized by experimental system and data modality to facilitate reproducibility and clarity.

---

## Repository structure
```
├── Scripts_multome_mm
├── Scripts_human_data
├── logs
└── README.md
```

## `Scripts_multiome_mm/`

Mouse endothelial cell **multiome single nucleus (GEX + ATAC)** analysis pipeline.

This folder contains sequentially numbered scripts covering:

- Multiome preprocessing and Seurat v5–based integration  
- Cell cycle scoring and endothelial subtype annotation  
- ATAC-seq peak calling and filtered peakset generation  
- Differential gene expression and chromatin accessibility analyses (edgeR)  
- Enhancer–gene linkage analysis  
- Pathway integration of RNA and ATAC data  
- Transcription factor motif enrichment (HOMER)  
- Heatmap and figure generation  

Scripts are numbered in the order they were executed to reproduce the analysis.

---

## `Scripts_human_data/`

Human endothelial cell **bulk RNA-seq, ATAC-seq, and ChIP-seq** analyses.

This folder includes scripts for:

- ATAC-seq peak merging and coverage quantification  
- Differential accessibility and expression analyses  
- Batch correction  
- ChIP-seq signal integration (e.g. JUNB, H3K27ac, H3K4me1)  
- Transcription factor motif enrichment  
- Human–mouse liftover and orthologous region processing  
- Cross-species correlation analyses between chromatin accessibility and gene expression  

---

## `logs/`

Log files generated during pipeline execution.

---

## Software and dependencies

Analyses were performed using a combination of R and command-line tools. Key dependencies include:

- R (≥ 4.2)  
- Seurat / Signac (v5)  
- edgeR  
- ComplexHeatmap  
- HOMER  
- bedtools  
- MACS2  

Individual scripts may require additional packages; see script headers for details.

---

## Notes

- Paths in the scripts assume a Linux-based environment.
- R scripts are meant to be run from the script folder.
- Bash scripts are meant to be run from main folder.
- Scripts were developed for internal reproducibility.  
- Raw and processed sequencing data are available through GEO (GSE313328).

---
