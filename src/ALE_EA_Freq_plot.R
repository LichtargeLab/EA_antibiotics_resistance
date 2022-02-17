# Author: Chen Wang, Dec 2021
# Graph ALE EA results

library(tidyverse)
library(scales)
source("src/functions.R")


cip.select <- tibble(GENE_name = c("gyrA", "parC", "parE", "ompF", "hemX", "udk", 
                                  "rstB", "yrbG", "sbcC", "mutL", "pdxY", "rob"),
                    class = c("both", "both", "both", "both", "control", "both",
                              "both", "Freq", "EA", "EA", "Freq", "Freq"))

# Graph cip EA-KS vs Freq
cip.graph.df <- openxlsx::read.xlsx("output/ALE_EA_output.xlsx", sheet = "cipro") 
cip.highlight.df <- cip.graph.df %>%
  left_join(cip.select) %>%
  filter(!is.na(class))
ggplot() +
  geom_point(data = filter(cip.graph.df, !GENE_name %in% cip.select$GENE_name), aes(x = Freq_combine, y = EAKS_combine), fill = "gray60", color = "gray60", pch=21, alpha = 0.5) +
  geom_point(data = cip.highlight.df, aes(x = Freq_combine, y = EAKS_combine, fill = class), color = "black",pch=21, size = 2) +
  ggrepel::geom_text_repel(data = cip.highlight.df, aes(x = Freq_combine, y = EAKS_combine, label = GENE_name), nudge_x = -0.28, nudge_y = 0, fontface = "bold", size = 4) +
  geom_function(fun = function(x) x, color = "red", linetype = 2) +
  scale_fill_manual(values = c("#FEF851", "gray10", "#418F21", "#EB58F9")) +
  guides(fill=guide_legend(title="")) +
  scale_x_continuous(trans=reverselog_trans(10)) +
  scale_y_continuous(trans=reverselog_trans(10)) +
  labs(x = "Frequency rank", y = "EA-KS rank") +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.title = element_text(face="bold", size = 15),
    axis.text = element_text(face="bold", size = 15),
    legend.text = element_text(face="bold", size = 10)
  ) 
ggsave("plot/ALE_EAvsFreq/ALE_cipro_KSvsFreq.jpeg", height = 1, width = 1.5, units = "in", scale = 4)


# Colistin
# Graph EA-KS vs Freq
col.select <- tibble(GENE_name = c("basS", "basR", "lapB", "waaQ", "lpxD", "ispB", "ybjX", "asmA",
                                  "ynjC", "yddW", "ddlB", "mntP", "srlA", "osmE"),
                    class = c("both", "both", "both", "both", "both", "both", "both", "both",
                              "EA", "EA", "EA", "Freq", "Freq", "Freq"))

col.graph.df <- openxlsx::read.xlsx("output/ALE_EA_output.xlsx", sheet = "colistin") 
col.highlight.df <- col.graph.df %>%
  left_join(col.select) %>%
  filter(!is.na(class))
ggplot() +
  geom_point(data = filter(col.graph.df, !GENE_name %in% col.select$GENE_name), aes(x = Freq_combine, y = EAKS_combine), fill = "gray60", color = "gray60", pch=21, alpha = 0.5) +
  geom_point(data = col.highlight.df, aes(x = Freq_combine, y = EAKS_combine, fill = class), color = "black",pch=21, size = 2) +
  ggrepel::geom_text_repel(data = col.highlight.df, aes(x = Freq_combine, y = EAKS_combine, label = GENE_name), nudge_x = -0.28, nudge_y = 0, fontface = "bold", size = 4) +
  geom_function(fun = function(x) x, color = "red", linetype = 2) +
  scale_fill_manual(values = c("#FEF851", "#418F21", "#EB58F9")) +
  guides(fill=guide_legend(title="")) +
  scale_x_continuous(trans=reverselog_trans(10)) +
  scale_y_continuous(trans=reverselog_trans(10)) +
  labs(x = "Frequency rank", y = "EA-KS rank") +
  coord_fixed(xlim = c(2300,1), ylim = c(2300,1)) +
  theme_classic() +
  theme(
    axis.title = element_text(face="bold", size = 15),
    axis.text = element_text(face="bold", size = 15),
    legend.text = element_text(face="bold", size = 10)
  )  
ggsave("plot/ALE_EAvsFreq/ALE_colistin_KSvsFreq.jpeg", height = 1, width = 1.5, units = "in", scale = 4)



## Graph EA vs Freq, labeling known genes

cip.known <- tibble(GENE_name = c("gyrA", "gyrB", "parC", "parE",
                                  "acrR", "marR", "soxR", "ompF"),
                    class = "known")

col.known <- tibble(GENE_name = c("basS", "basR"),
                    class = "known")

GraphRankings(cip.graph.df, cip.known, "Freq_combine", "EAKS_combine", 
              title = NULL, ylab = "EA-KS rank") +
  theme(legend.position = "none")
ggsave("plot/ALE_EAvsFreq/ALE_cipro_EAKSvsFreq_known.jpeg", height = 1, width = 1.5, units = "in", scale = 4)

GraphRankings(cip.graph.df, cip.known, "Freq_combine", "EAsum_combine", 
              title = NULL, ylab = "EA-sum rank") +
  theme(legend.position = "none")
ggsave("plot/ALE_EAvsFreq/ALE_cipro_EAsumvsFreq_known.jpeg", height = 1, width = 1.5, units = "in", scale = 4)

GraphRankings(col.graph.df, col.known, "Freq_combine", "EAKS_combine", 
              title = NULL, ylab = "EA-KS rank") +
  theme(legend.position = "none")
ggsave("plot/ALE_EAvsFreq/ALE_colistin_EAKSvsFreq_known.jpeg", height = 1, width = 1.5, units = "in", scale = 4)

GraphRankings(col.graph.df, col.known, "Freq_combine", "EAsum_combine", 
              title = NULL, ylab = "EA-sum rank") +
  theme(legend.position = "none")
ggsave("plot/ALE_EAvsFreq/ALE_colistin_EAsumvsFreq_known.jpeg", height = 1, width = 1.5, units = "in", scale = 4)




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


cip_by_condition <- tibble(x_var = c("Freq.rank"), y_var = c("EAKS.rank", "EAsum.rank"),
                           x_lab = c("Frequency rank"), y_lab = c("EA_KS rank", "EA_sum rank")) %>%
  crossing(., tibble(mut = c("MT", "WM", "WT"),
                     mut_full = c("mutator", "WT+mutagen", "WT"))) %>%
  mutate(x_var = paste0(x_var, "_", mut),
         y_var = paste0(y_var, "_", mut)) %>%
  mutate(fig = pmap(list(x_var, y_var, x_lab, y_lab, mut_full),
                    ~GraphRankings(cip.graph.df, cip.genes, ..1, ..2, ..5, ..3,..4) +
                      theme(legend.position = "none"))) %>%
  arrange((y_lab), desc(mut))

cowplot::plot_grid(plotlist = cip_by_condition$fig, align = "hv", nrow = 2)


col_by_condition <- tibble(x_var = c("Freq.rank"), y_var = c("EAKS.rank", "EAsum.rank"),
                           x_lab = c("Frequency rank"), y_lab = c("EA_KS rank", "EA_sum rank")) %>%
  crossing(., tibble(mut = c("MT", "WM", "WT"),
                     mut_full = c("mutator", "WT+mutagen", "WT"))) %>%
  mutate(x_var = paste0(x_var, "_", mut),
         y_var = paste0(y_var, "_", mut)) %>%
  mutate(fig = pmap(list(x_var, y_var, x_lab, y_lab, mut_full),
                    ~GraphRankings(col.graph.df, col.genes, ..1, ..2, ..5, ..3,..4) +
                      theme(legend.position = "none"))) %>%
  arrange((y_lab), desc(mut))

cowplot::plot_grid(plotlist = col_by_condition$fig, align = "hv", nrow = 2)
