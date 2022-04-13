# Author: Chen Wang, Dec 2021
# Scatter column plot for ALE

library(tidyverse)
library(scales)
library(ggbeeswarm)
source("src/functions.R")


cip.graph.df <- openxlsx::read.xlsx("output/ALE_EA_output.xlsx", sheet = "cipro") 
col.graph.df <- openxlsx::read.xlsx("output/ALE_EA_output.xlsx", sheet = "colistin") 

cip.known <- tibble(GENE_name = c("gyrA", "gyrB", "parC", "parE",
                                  "acrR", "marR", "soxR", "ompF"),
                    class = "known")

col.known <- tibble(GENE_name = c("basS", "basR"),
                    class = "known")

cip.test <- cip.graph.df %>%
  select(GENE_name, starts_with("EAKS.rank_"), starts_with("EAsum.rank_")) %>%
  pivot_longer(cols = 2:7, names_to = "type", values_to = "rank") %>%
  separate(col = type, into = c("method", "mut"), sep = "_") %>%
  mutate(method = if_else(method == "EAsum.rank", "EA_sum", "EA_KS")) %>%
  mutate(mut = factor(mut, c("WT", "WM", "MT")),
         mut = fct_recode(mut,
                          `WT+mutagen` = "WM",
                          mutator = "MT")) %>%
  left_join(cip.known) %>%
  mutate(label = ifelse(class == "known", GENE_name, NA)) %>%
  replace_na(list(class = "other")) %>%
  arrange(desc(class)) 


ggplot(data = cip.test) +
  geom_quasirandom(aes(x = " ", y = rank, fill = class, alpha = class, size = class, color = class),
                   shape = 21, show.legend = FALSE,
                   position = position_quasirandom(bandwidth = 0.8)) +
  ggrepel::geom_text_repel(aes(x = " ", y = rank, label = label),
                           position = position_quasirandom(bandwidth = 0.8), fontface = "italic") +
  scale_y_continuous(trans=reverselog_trans(10), minor_breaks = NULL) +
  scale_color_manual(values = c("black", "grey70")) +
  scale_fill_manual(values = c("dark orange", "white")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_size_manual(values = c(2, 1.5)) +
  facet_grid(mut~method, scales = "free") +
  xlab("") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11), # set font for label bars
        axis.text = element_text(size = 11), # set font for axis numbers
        axis.title = element_text(size = 11), # set font for axis titles
        title = element_text(size = 11))
ggsave("plot/ALE_EAvsFreq/ALE_EA_scatter_column_cipro.pdf", width = 4, height = 6, units = "in")


col.test <- col.graph.df %>%
  select(GENE_name, starts_with("EAKS.rank_"), starts_with("EAsum.rank_")) %>%
  pivot_longer(cols = 2:7, names_to = "type", values_to = "rank") %>%
  separate(col = type, into = c("method", "mut"), sep = "_") %>%
  mutate(method = if_else(method == "EAsum.rank", "EA_sum", "EA_KS")) %>%
  mutate(mut = factor(mut, c("WT", "WM", "MT")),
         mut = fct_recode(mut,
                          `WT+mutagen` = "WM",
                          mutator = "MT")) %>%
  left_join(col.known) %>%
  mutate(label = ifelse(class == "known", GENE_name, NA)) %>%
  replace_na(list(class = "other")) %>%
  arrange(desc(class)) 

ggplot(data = col.test) +
  geom_quasirandom(aes(x = " ", y = rank, fill = class, alpha = class, size = class, color = class),
                   shape = 21, show.legend = FALSE,
                   position = position_quasirandom(bandwidth = 0.8)) +
  ggrepel::geom_text_repel(aes(x = " ", y = rank, label = label),
                           position = position_quasirandom(bandwidth = 0.8), fontface = "italic") +
  scale_y_continuous(trans=reverselog_trans(10), minor_breaks = NULL) +
  scale_color_manual(values = c("black", "grey70")) +
  scale_fill_manual(values = c("dark orange", "white")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_size_manual(values = c(2, 1.5)) +
  facet_grid(mut~method, scales = "free") +
  xlab("") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11), # set font for label bars
        axis.text = element_text(size = 11), # set font for axis numbers
        axis.title = element_text(size = 11), # set font for axis titles
        title = element_text(size = 11))
ggsave("plot/ALE_EAvsFreq/ALE_EA_scatter_column_colistin.pdf", width = 4, height = 6, units = "in")
