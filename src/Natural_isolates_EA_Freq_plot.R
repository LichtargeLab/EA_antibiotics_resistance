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

cip_EA <- openxlsx::read.xlsx("output/Natural_isolates_EA_output.xlsx", sheet = "cipro") 


col_EA <- openxlsx::read.xlsx("output/Natural_isolates_EA_output.xlsx", sheet = "colistin") 





GraphRankings(cip_EA, highlight.genes = cip.genes, xvar = "Freq.rank", yvar = "EAKS.rank", 
              title = NULL,
              xlab = "Frequency rank", ylab = "EA_KS rank") +
  theme(legend.position = "none")
ggsave("plot/Natural_isolates_EAvsFreq/Natural_isolates_cipro_EAKSvsFreq.pdf", height = 1, width = 1.2, units = "in", scale = 4)

GraphRankings(cip_EA, highlight.genes = cip.genes, xvar = "Freq.rank", yvar = "EAsum.rank", 
              title = NULL,
              xlab = "Frequency rank", ylab = "EA_sum rank") +
  theme(legend.position = "none")
ggsave("plot/Natural_isolates_EAvsFreq/Natural_isolates_cipro_EAsumvsFreq.pdf", height = 1, width = 1.2, units = "in", scale = 4)



GraphRankings(col_EA, highlight.genes = col.genes, xvar = "Freq.rank", yvar = "EAKS.rank", 
              title = NULL,
              xlab = "Frequency rank", ylab = "EA_KS rank") +
  theme(legend.position = "none")
ggsave("plot/Natural_isolates_EAvsFreq/Natural_isolates_colistin_EAKSvsFreq.pdf", height = 1, width = 1.2, units = "in", scale = 4)

GraphRankings(col_EA, highlight.genes = col.genes, xvar = "Freq.rank", yvar = "EAsum.rank", 
              title = NULL,
              xlab = "Frequency rank", ylab = "EA_sum rank") +
  theme(legend.position = "none")
ggsave("plot/Natural_isolates_EAvsFreq/Natural_isolates_colistin_EAsumvsFreq.pdf", height = 1, width = 1.2, units = "in", scale = 4)


legend <- cowplot::get_legend(GraphRankings(col_EA, highlight.genes = col.genes, xvar = "Freq.rank", yvar = "EAsum.rank", 
                                            title = NULL,
                                            xlab = "Frequency rank", ylab = "EA_sum rank"))

pdf("plot/Natural_isolates_EAvsFreq/legend.pdf", height = 0.8, width = 2)
grid::grid.newpage()
grid::grid.draw(legend)
dev.off()
