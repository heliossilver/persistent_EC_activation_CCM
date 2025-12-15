
pacman::p_load("purrr", "dplyr", "tidyr", "glue", "stringr")

ortholog_dir <- "../results/human_mouse_integration/data"
dir.create(ortholog_dir, recursive = TRUE, showWarnings = FALSE)


files <- list.files("../results/human_mouse_integration/liftover/", pattern = "intersect.tsv")

df <- tibble(file    = files,
             sample_names  = str_extract(files, "^[^_]+_[^_]+_[^_]+"),  # e.g., D3_vs_Mm
             names    = str_extract(files, "^[^_]+"))  %>%              # e.g., D3
  distinct(sample_names, names)

pwalk(df, function(sample_names, names) {
  # A: human->mouse (BtoA_vs_A)
  human_to_mouse <- read.table(
    glue("../results/human_mouse_integration/liftover/{sample_names}_AtoB_vs_B.intersect.tsv"),
    sep = "\t", header = FALSE, quote = "", comment.char = "") %>%
    arrange(desc(V9)) %>%
    distinct(V8, .keep_all = TRUE) %>%
    select(V4, V8) %>%
    dplyr::rename(peaks_hu = V4, peaks_ms = V8)
  
  # B: mouse->human (AtoB_vs_B)
  mouse_to_human <- read.table(
    glue("../results/human_mouse_integration/liftover/{sample_names}_BtoA_vs_A.intersect.tsv"),
    sep = "\t", header = FALSE, quote = "", comment.char = "") %>%
    arrange(desc(V9)) %>%
    distinct(V8, .keep_all = TRUE) %>%
    select(V4, V8) %>%
    dplyr::rename(peaks_ms = V4, peaks_hu = V8)
  
  combined_unique <- bind_rows(human_to_mouse, mouse_to_human) %>%
    distinct(peaks_hu, peaks_ms, .keep_all = TRUE) %>%
    select(peaks_hu, peaks_ms)
  
  write.csv(combined_unique,
            glue("{ortholog_dir}/{names}_hu_ms_ortholog_peaks.csv"),
            row.names = FALSE)
})
