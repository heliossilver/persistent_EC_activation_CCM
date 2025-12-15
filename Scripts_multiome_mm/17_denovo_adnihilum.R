
pacman::p_load("dplyr", "tibble", "tidyr", "readr", "purrr",
                "stringr", "rtracklayer", "GenomicRanges", 
                 "glue", "ggplot2", "forcats")



### De Novo peaks----
ord <- c(1, 3, 2, 4, 5)

labs <- c("Artery KO", "Capillary artery KO", "Capillary KO", "Capillary vein KO", "Large vein KO") # Build summary table 
peaks_df_l <- list.files("../results/ATAC/edgeR/dar_diffacc_tables/", pattern = "_vs_Flox_diffacc.csv", full.names = TRUE) %>% 
  map(~ vroom::vroom(.x, show_col_types = FALSE) %>% 
        dplyr::filter(diffaccessible == "Up") %>% 
        dplyr::select(peaks) %>% 
        mutate(peak_label = peaks) %>% 
        separate(peaks, into = c("chr", "start", "end"), sep = "-", convert = TRUE) %>%
        makeGRangesFromDataFrame(keep.extra.columns = T)) %>% 
  magrittr::extract(ord) %>% set_names(labs) 

# 2) Flox peaksets as GRanges per subtype (explicitly target the Flox files) 

flox_files <- list.files("../results/ATAC/merged_peaks/", pattern = "_Flox\\.filterNfixed\\.peakset$", full.names = TRUE) 
gr_flox_df_l <- flox_files %>% 
  map(~ vroom::vroom(.x, show_col_types = FALSE) %>% 
        filter(spm >= 2) %>% 
        makeGRangesFromDataFrame(keep.extra.columns = T)) %>% 
  set_names(labs) 


gr_denovo <- map2( peaks_df_l, gr_flox_df_l,  ~ { 
  subsetByOverlaps(.x, .y, invert = TRUE) %>% 
    GenomicRanges::reduce()  }) 


dir.create("../results/RNA_ATAC_integration/denovo_cCRE_KO_per_cell", showWarnings = FALSE, recursive = TRUE) 

walk2(gr_denovo, names(gr_denovo), ~ { 
  export(.x, glue::glue("../results/RNA_ATAC_integration/denovo_cCRE_KO_per_cell/{gsub(' ', '_', .y)}_denovo.bed")) }) 

dar_df_l <- list.files("../results/ATAC/edgeR/dar_diffacc_tables/", pattern = "_vs_Flox_diffacc.csv", full.names = TRUE) %>% 
  map(~ vroom::vroom(.x, show_col_types = FALSE) %>% 
        dplyr::filter(diffaccessible == "Up") %>% 
        dplyr::select(peaks, starts_with("logFC"), starts_with("FDR")) %>% 
  dplyr::rename_with(~ paste0(.x, "_peak"), starts_with("logFC")) %>%
  dplyr::rename_with(~ paste0(.x, "_peak"), starts_with("FDR"))) %>% 
  magrittr::extract(ord) %>% set_names(labs) 

deg_df_l <- list.files("../results/RNA/edgeR/deg_difexp_tables/", pattern = "_vs_Flox_diffexp.csv", full.names = TRUE) %>% 
  map(~ vroom::vroom(.x, show_col_types = FALSE) %>% 
        dplyr::select(gene, starts_with("logFC"), starts_with("FDR")) %>% 
        dplyr::rename_with(~ paste0(.x, "_gene"), starts_with("logFC")) %>%
        dplyr::rename_with(~ paste0(.x, "_gene"), starts_with("FDR"))) %>% 
  magrittr::extract(ord) %>% set_names(labs) 

links_df <- read.csv("../results/ATAC/data/endo_linkpeaks.csv") %>% 
  dplyr::select(peaks, gene) 

gr_denovo_df <- map2(gr_denovo, names(dar_df_l), ~{ 
  deg_df <- deg_df_l[[.y]]
  as.data.frame(.x) %>% unite(peaks, seqnames, start, end, sep = "-") %>% 
    select(-c(width, strand)) %>% 
    inner_join(dar_df_l[[.y]], by = "peaks") %>% 
    left_join(links_df, by = "peaks") %>% 
    mutate(putative_enhancer = factor(if_else(!is.na(gene), "Yes", "No"), 
                                      levels = c("Yes","No")), 
           .logFC = pick(starts_with("logFC"))[[1]]) %>%
    arrange(putative_enhancer, 
            desc(.logFC)) %>% dplyr::select(-.logFC) %>% 
    left_join(deg_df, by = "gene")
  }) 

walk2(gr_denovo_df, names(gr_denovo), ~ { 
  write_csv(.x, glue::glue("../results/RNA_ATAC_integration/denovo_cCRE_KO_per_cell/{gsub(' ', '_', .y)}_denovo.csv")) })

# a <- purrr::reduce(gr_denovo_df, full_join, by = c("peaks", "putative_enhancer", "gene")) %>% 
#   distinct(peaks, gene, .keep_all = TRUE)

### Ad nihilum peaks----

peaks_down_df_l <- list.files("../results/ATAC/edgeR/dar_diffacc_tables/", pattern = "_vs_Flox_diffacc.csv", full.names = TRUE) %>% 
  map(~ vroom::vroom(.x, show_col_types = FALSE) %>% 
        dplyr::filter(diffaccessible == "Down") %>% 
        dplyr::select(peaks) %>% 
        mutate(peak_label = peaks) %>% 
        separate(peaks, into = c("chr", "start", "end"), sep = "-", convert = TRUE) %>%
        makeGRangesFromDataFrame(keep.extra.columns = T)) %>% 
  magrittr::extract(ord) %>% 
  set_names(labs) 


ko_files <- list.files("../results/ATAC/merged_peaks/", pattern = "_KO\\.filterNfixed\\.peakset$", full.names = TRUE) 
gr_ko_df_l <- ko_files %>% 
  map(~ vroom::vroom(.x, show_col_types = FALSE) %>% 
        filter(spm >= 2) %>% 
        makeGRangesFromDataFrame(keep.extra.columns = T)) %>% 
  set_names(labs) 


gr_adnihilum <- map2( peaks_down_df_l, gr_ko_df_l,  ~ { 
  subsetByOverlaps(.x, .y, invert = TRUE) %>% 
    GenomicRanges::reduce()  }) 


dir.create("../results/RNA_ATAC_integration/adnihilum_cCRE_KO_per_cell", showWarnings = FALSE, recursive = TRUE) 

walk2(gr_adnihilum, names(gr_adnihilum), ~ { 
  export(.x, glue::glue("../results/RNA_ATAC_integration/adnihilum_cCRE_KO_per_cell/{.y}_adnihilum.bed")) }) 

dar_down_df_l <- list.files("../results/ATAC/edgeR/dar_diffacc_tables/", pattern = "_vs_Flox_diffacc.csv", full.names = TRUE) %>% 
  map(~ vroom::vroom(.x, show_col_types = FALSE) %>% 
        dplyr::filter(diffaccessible == "Down") %>% 
        dplyr::select(peaks, starts_with("logFC"), starts_with("FDR"))) %>% 
  magrittr::extract(ord) %>% set_names(labs) 


links_df <- read.csv("../results/ATAC/data/endo_linkpeaks.csv") %>% 
  dplyr::select(peaks, gene) 

gr_adnihilum_df <- map2(gr_adnihilum, names(dar_down_df_l), ~{ 
  deg_df <- deg_df_l[[.y]]
  as.data.frame(.x) %>% unite(peaks, seqnames, start, end, sep = "-") %>% 
    select(-c(width, strand)) %>% 
    inner_join(dar_down_df_l[[.y]], by = "peaks") %>% 
    left_join(links_df, by = "peaks") %>%
    mutate(putative_enhancer = factor(if_else(!is.na(gene), "Yes", "No"),
                                      levels = c("Yes","No")),
           .logFC = pick(starts_with("logFC"))[[1]]) %>% # first column that starts with "logFC"
    arrange(putative_enhancer, .logFC) %>% 
    dplyr::select(-.logFC )%>% 
    left_join(deg_df, by = "gene")
  })

walk2(gr_adnihilum_df, names(gr_adnihilum_df), ~ { 
  write_csv(.x, glue::glue("../results/RNA_ATAC_integration/adnihilum_cCRE_KO_per_cell/{.y}_adnihilum.csv")) })

###summary de novo----

n_denovo <- map_dbl(gr_denovo_df, ~{
  .x %>% distinct(peaks) %>% nrow()
})

n_enhancers <- map_dbl(gr_denovo_df, ~{
  .x %>% filter(putative_enhancer == "Yes") %>% distinct(peaks) %>% nrow()
})


denovo_summ <- tibble(cell_type = fct_rev(labs),
                      n_denovo = n_denovo,
                      n_enhancers = n_enhancers)

denovo_long <- denovo_summ %>%
  mutate(n_non_enhancers = n_denovo - n_enhancers,  
         cell_type = factor(cell_type, levels = rev(cell_type))) %>%
  select(cell_type,
         `Putative enhancers` = n_enhancers,
         `Non-putative enhancers` = n_non_enhancers) %>%
  pivot_longer(cols = c(`Putative enhancers`, `Non-putative enhancers`),   
               names_to = "type",
               values_to = "n") %>% 
  mutate(type = factor(type, levels = c("Putative enhancers", "Non-putative enhancers")))

###summary ad nihilum----
n_adnihilum <- map_dbl(gr_adnihilum_df, ~{
  .x %>% distinct(peaks) %>% nrow()
})

n_enhancers_adnihilum <- map_dbl(gr_adnihilum_df, ~{
  .x %>% filter(putative_enhancer == "Yes") %>% distinct(peaks) %>% nrow()
})

adnihilum_summ <- tibble(cell_type = fct_rev(labs),
                      n_adnihilum = n_adnihilum,
                      n_enhancers = n_enhancers_adnihilum)

adnihilum_long <- adnihilum_summ %>%
  mutate(n_non_enhancers = n_adnihilum - n_enhancers,  
         cell_type = factor(cell_type, levels = rev(cell_type))) %>%
  select(cell_type,
         `Putative enhancers` = n_enhancers,
         `Non-putative enhancers` = n_non_enhancers) %>%
  pivot_longer(cols = c(`Putative enhancers`, `Non-putative enhancers`),   
               names_to = "type",
               values_to = "n") %>% 
  mutate(type = factor(type, levels = c("Putative enhancers", "Non-putative enhancers")))


### plots ----
denovo_col <- RColorBrewer::brewer.pal(12, "Set3")[3]
adnihi_col <- RColorBrewer::brewer.pal(12, "Set3")[8]

EC_colors <- rev(c("#247567", "#E2705B", "#6C7A99", "#B4698B", "#D8A21F"))

denovo_plot <- ggplot(denovo_summ, aes(x = cell_type, y = n_denovo)) +
  geom_col(colour = denovo_col, fill = denovo_col) +
  geom_text(data = denovo_summ,
            aes(x = cell_type, y = n_denovo, label = n_denovo),
            inherit.aes = FALSE,
            hjust = -0.2,           # because of coord_flip()
            size = 3.5,
            fontface = "bold") + 
  labs(x = NULL,
       y = "Number of de novo cCREs",
       fill = NULL,
       title = "<b>De novo cCREs per BEC subtype</b><br><b>(</b><i><b>Pdcd10</b></i><sup><b>BECKO</b></sup> vs <i><b>Pdcd10</b></i><sup><b>fl/fl</b></sup><b>)</b>") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(face = "bold", margin = margin(t = 10, unit = "pt")),
        axis.text.y = element_text(size = 8, hjust = 1, face = "bold", colour = EC_colors),
        plot.title = ggtext::element_markdown(hjust = 0.5, size = 10),
        legend.text = element_text(size = 8, face = "bold")) + coord_flip(clip = "off") +
  scale_y_continuous(breaks =  seq(0, 14000, by = 4000), labels = function(x) abs(x), limits = c(0, 17000))


adnihilum_plot <- ggplot(adnihilum_summ, aes(x = cell_type, y = n_adnihilum)) +
  geom_col(colour = "#F4C2C2", fill = "#F4C2C2") +
  geom_text(data = adnihilum_summ,
            aes(x = cell_type, y = n_adnihilum, label = n_adnihilum),
            inherit.aes = FALSE,
            hjust = -0.2,           # because of coord_flip()
            size = 3.5,
            fontface = "bold") + 
  labs(x = NULL,
       y = "Number of ad nihilum cCREs",
       fill = NULL,
       title = "<b>Ad nihilum cCREs per BEC subtype</b><br><b>(</b><i><b>Pdcd10</b></i><sup><b>BECKO</b></sup> vs <i><b>Pdcd10</b></i><sup><b>fl/fl</b></sup><b>)</b>") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text( hjust = 0.5, size = 9),
        axis.title.x = element_text(face = "bold", margin = margin(t = 10, unit = "pt")),
        axis.text.y = element_text(size = 8, hjust = 1, face = "bold", colour = EC_colors),
        plot.title = ggtext::element_markdown(hjust = 0.5, size = 10),
        legend.text = element_text(size = 8, face = "bold")) + coord_flip(clip = "off") +
  scale_y_continuous(breaks =  seq(0, 4000, by = 1500), labels = function(x) abs(x), limits = c(0, 4000))

ggsave("../results/RNA_ATAC_integration/graphs/cCRE_KO_de_novo.png", denovo_plot, width = 3.5, height = 2.5, units = "in", dpi = 300)
ggsave("../results/RNA_ATAC_integration/graphs/cCRE_KO_ad_nihilum.png", adnihilum_plot, width = 3.5, height = 2.5, units = "in", dpi = 300)



fraction_denovo <- ggplot(denovo_long, aes(x = cell_type, y = n, fill = type )) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c( "Non-putative enhancers" = denovo_col,
                    "Putative enhancers" = RColorBrewer::brewer.pal(12, "Set3")[10]),
                    guide = guide_legend(reverse = TRUE)) +
  geom_text(data = denovo_summ %>%
              mutate(non_enhancer_frac = 1 - (n_enhancers / n_denovo),
                     position = non_enhancer_frac + ((1 - non_enhancer_frac)/2)),
            aes(x = cell_type, y = position, label = n_enhancers),
            inherit.aes = FALSE,
            size = 2.5,
            colour = "white",
            fontface = "bold") +  
  geom_text(data = denovo_summ %>%
              mutate(non_enhancer_frac = 1 - (n_enhancers / n_denovo),
                     position = non_enhancer_frac/2),
            aes(x = cell_type, y = position, label = n_denovo - n_enhancers),
            inherit.aes = FALSE,
            size = 2.5,
            colour = "white",
            fontface = "bold") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 8, hjust = 1, face = "bold", color = EC_colors),
        legend.position = "top",
        legend.direction = "horizontal") + 
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL,
       y = "Fraction of de novo cCREs",
       fill = NULL)  + coord_flip() 

# ggplot_build(fraction_denovo)$layout$panel_params[[1]]$y$get_labels()

fraction_ad_nihilum <- ggplot(adnihilum_long, aes(x = cell_type, y = n, fill = type )) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("Putative enhancers" = RColorBrewer::brewer.pal(8, "Set1")[8],        # red / accent
                               "Non-putative enhancers" = "#F4C2C2"),
                    guide = guide_legend(reverse = TRUE)) +
  geom_text(data = adnihilum_summ %>%
              mutate(non_enhancer_frac = 1 - (n_enhancers / n_adnihilum),
                     position = non_enhancer_frac + ((1 - non_enhancer_frac)/2)),
            aes(x = cell_type, y = position, label = n_enhancers),
            inherit.aes = FALSE,
            size = 2.5,
            colour = "white", 
            fontface = "bold") +  
  geom_text(data = adnihilum_summ %>%
              mutate(non_enhancer_frac = 1 - (n_enhancers / n_adnihilum),
                     position = non_enhancer_frac/2 ),
            aes(x = cell_type, y = position, label = n_adnihilum - n_enhancers),
            inherit.aes = FALSE,
            size = 2.5,
            colour = "white",
            hjust = -0.1,  
            fontface = "bold")  + 
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL,
       y = "Fraction of ad nihilum cCREs",
       fill = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        axis.text.y = element_text(size = 8, hjust = 1, face = "bold", colour = EC_colors)) + coord_flip() 


ggsave("../results/RNA_ATAC_integration/graphs/cCRE_KO_de_novo_fractions.png", fraction_denovo, width = 4, height = 2.5, units = "in", dpi = 300)
ggsave("../results/RNA_ATAC_integration/graphs/cCRE_KO_ad_nihilum_fractions.png", fraction_ad_nihilum, width = 4, height = 2.5, units = "in", dpi = 300)


# === Log R session info for reproducibility ===
sink("../logs/sessioninfo_17_denovo_adnihilum.txt")
sessionInfo()
sink()

