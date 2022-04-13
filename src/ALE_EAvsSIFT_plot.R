# Author: Chen Wang, Dec 2021

library(tidyverse)
source("src/functions.R")


cip.genes <- tibble(GENE_name = c("gyrA", "gyrB", "parC", "parE",
                                  "acrR", "marR", "soxR", "ompF"),
                    class = "known") %>%
  bind_rows(tibble(GENE_name = c("rob"),
                   class = "newly_identified"))

col.genes <- tibble(GENE_name = c("basS", "basR"),
                    class = "known") %>%
  bind_rows(tibble(GENE_name = c("asmA", "lapB", "ispB", "waaQ",
                                 "ybjX", "ynjC", "osmE"),
                   class = "newly_identified"))

cip_EA <- openxlsx::read.xlsx("output/ALE_EA_output.xlsx", sheet = "cipro") %>%
  select(b_num, GENE_name, contains("EAKS"), contains("EAsum"))

cip_SIFT <- openxlsx::read.xlsx("output/ALE_SIFT_output.xlsx", sheet = "cipro") %>%
  select(b_num, GENE_name, contains("SIFT_adj_KS"), contains("SIFT_adj_sum"))

col_EA <- openxlsx::read.xlsx("output/ALE_EA_output.xlsx", sheet = "colistin") %>%
  select(b_num, GENE_name, contains("EAKS"), contains("EAsum"))

col_SIFT <- openxlsx::read.xlsx("output/ALE_SIFT_output.xlsx", sheet = "colistin") %>%
  select(b_num, GENE_name, contains("SIFT_adj_KS"), contains("SIFT_adj_sum"))


cip <- full_join(cip_EA, cip_SIFT)
col <- full_join(col_EA, col_SIFT)

GraphRankings(cip, highlight.genes = cip.genes, xvar = "SIFT_adj_KS_combine", yvar = "EAKS_combine", 
              title = NULL,
              xlab = "SIFT_adj_KS rank", ylab = "EA_KS rank")
ggsave("plot/ALE_EAvsSIFT/ALE_cipro_EAKSvsSIFTKS.pdf", height = 1, width = 1.5, units = "in", scale = 4)


GraphRankings(cip, highlight.genes = cip.genes, xvar = "SIFT_adj_sum_combine", yvar = "EAsum_combine", 
              title = NULL, 
              xlab = "SIFT_adj_sum rank", ylab = "EA_sum rank")
ggsave("plot/ALE_EAvsSIFT/ALE_cipro_EAsumvsSIFTsum.pdf", height = 1, width = 1.5, units = "in", scale = 4)



GraphRankings(col, highlight.genes = col.genes, xvar = "SIFT_adj_KS_combine", yvar = "EAKS_combine", 
              title = NULL, 
              xlab = "SIFT_adj_KS rank", ylab = "EA_KS rank")
ggsave("plot/ALE_EAvsSIFT/ALE_colistin_EAKSvsSIFTKS.pdf", height = 1, width = 1.5, units = "in", scale = 4)


GraphRankings(col, highlight.genes = col.genes, xvar = "SIFT_adj_sum_combine", yvar = "EAsum_combine", 
              title = NULL, 
              xlab = "SIFT_adj_sum rank", ylab = "EA_sum rank")
ggsave("plot/ALE_EAvsSIFT/ALE_colistin_EAsumvsSIFTsum.pdf", height = 1, width = 1.5, units = "in", scale = 4)








cip_by_condition <- tibble(x_var = c("SIFT_adj_KS.rank", "SIFT_adj_sum.rank"), y_var = c("EAKS.rank", "EAsum.rank"),
                           x_lab = c("SIFT_adj KS rank", "SIFT_adj sum rank"), y_lab = c("EA_KS rank", "EA_sum rank")) %>%
  crossing(., tibble(mut = c("MT", "WM", "WT"),
                     mut_full = c("mutator", "WT+mutagen", "WT"))) %>%
  mutate(x_var = paste0(x_var, "_", mut),
         y_var = paste0(y_var, "_", mut)) %>%
  mutate(fig = pmap(list(x_var, y_var, x_lab, y_lab, mut_full),
                    ~GraphRankings(cip, cip.genes, ..1, ..2, ..5, ..3,..4) +
                      theme(legend.position = "none"))) %>%
  arrange((y_lab), desc(mut))

cowplot::plot_grid(plotlist = cip_by_condition$fig, align = "hv", nrow = 2)


col_by_condition <- tibble(x_var = c("SIFT_adj_KS.rank", "SIFT_adj_sum.rank"), y_var = c("EAKS.rank", "EAsum.rank"),
                           x_lab = c("SIFT_adj KS rank", "SIFT_adj sum rank"), y_lab = c("EA_KS rank", "EA_sum rank")) %>%
  crossing(., tibble(mut = c("MT", "WM", "WT"),
                     mut_full = c("mutator", "WT+mutagen", "WT"))) %>%
  mutate(x_var = paste0(x_var, "_", mut),
         y_var = paste0(y_var, "_", mut)) %>%
  mutate(fig = pmap(list(x_var, y_var, x_lab, y_lab, mut_full),
                    ~GraphRankings(col, col.genes, ..1, ..2, ..5, ..3,..4) +
                      theme(legend.position = "none"))) %>%
  arrange((y_lab), desc(mut))
cowplot::plot_grid(plotlist = col_by_condition$fig, align = "hv", nrow = 2)
