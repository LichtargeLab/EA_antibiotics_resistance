# Author: Chen Wang, Dec 2021
# Graph natural isolates EA results

library(tidyverse)
library(scales)
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

cip_EA <- openxlsx::read.xlsx("output/Natural_isolates_EA_output.xlsx", sheet = "cipro") %>%
  select(-Freq.rank)

col_EA <- openxlsx::read.xlsx("output/Natural_isolates_EA_output.xlsx", sheet = "colistin") %>%
  select(-Freq.rank)

cip_SIFT <- openxlsx::read.xlsx("output/Natural_isolates_SIFT_output.xlsx", sheet = "cipro") %>%
  select(-Freq.rank)

col_SIFT <- openxlsx::read.xlsx("output/Natural_isolates_SIFT_output.xlsx", sheet = "colistin") %>%
  select(-Freq.rank)


cip <- full_join(cip_EA, cip_SIFT)
col <- full_join(col_EA, col_SIFT)


openxlsx::write.xlsx(list(cip, col),
                     "output/Natural_isolates_EAvsSIFT.xlsx", sheetName = c("cipro", "colistin"))


GraphRankings(cip, highlight.genes = cip.genes, xvar = "SIFT_adj_KS.rank", yvar = "EAKS.rank", 
              title = NULL,
              xlab = "SIFT_adj_KS rank", ylab = "EA_KS rank")
ggsave("plot/Natural_isolates_EAvsSIFT/Natural_isolates_cipro_EAKSvsSIFTKS.pdf", height = 1, width = 1.5, units = "in", scale = 4)

GraphRankings(cip, highlight.genes = cip.genes, xvar = "SIFT_adj_sum.rank", yvar = "EAsum.rank", 
              title = NULL,
              xlab = "SIFT_adj_sum rank", ylab = "EA_sum rank")
ggsave("plot/Natural_isolates_EAvsSIFT/Natural_isolates_cipro_EAsumvsSIFTsum.pdf", height = 1, width = 1.5, units = "in", scale = 4)



GraphRankings(col, highlight.genes = col.genes, xvar = "SIFT_adj_KS.rank", yvar = "EAKS.rank", 
              title = NULL,
              xlab = "SIFT_adj_KS rank", ylab = "EA_KS rank")
ggsave("plot/Natural_isolates_EAvsSIFT/Natural_isolates_colistin_EAKSvsSIFTKS.pdf", height = 1, width = 1.5, units = "in", scale = 4)

GraphRankings(col, highlight.genes = col.genes, xvar = "SIFT_adj_sum.rank", yvar = "EAsum.rank", 
              title = NULL,
              xlab = "SIFT_adj_sum rank", ylab = "EA_sum rank")
ggsave("plot/Natural_isolates_EAvsSIFT/Natural_isolates_colistin_EAsumvsSIFTsum.pdf", height = 1, width = 1.5, units = "in", scale = 4)
