#!/bin/bash

### After calling peaks with MACS2, summit paths were included in a table named folder_summit_paths.txt

Rscript iterative_overlap_peak_merging.R \
	-i folder_summit_paths.txt \
	-g hg38 \
	-d ~/mBEC_multiome/Human_D3/D3_ATAC/data/D3_merged_peaks/ \
	-o condition_treatment



