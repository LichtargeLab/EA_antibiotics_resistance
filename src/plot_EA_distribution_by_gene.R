# Author: Chen Wang, Dec 2021
# Plot the EA distribution for a given gene and compare against random
# mutation background.
library(tidyverse)
library(stickylabeller)

len_map <- read_csv("data/len_map.csv")

all_data <- read_csv("data/ALE_mutations_with_SIFT.csv",
                     col_types = "dccccddcccd")

all_data <- all_data %>%
  filter(!is.na(SIFT_adj))

random <- filter(all_data, mut == "random") %>%
  filter(ACTION > 0 & ACTION <= 100)

d0 <- filter(all_data, abx == "d0") %>%
  filter(ACTION > 0 & ACTION <= 100)

evolve <- filter(all_data, abx %in% c("CIP", "COL")) %>% 
  left_join(select(len_map, GN = GENE_name, b_num), by = "b_num") %>%
  anti_join(., d0, by = c("POS", "REF", "ALT", "mut")) %>%
  filter(ACTION > 0 & ACTION <= 100) %>%
  filter(AF > 0.1) %>%
  mutate(mut = factor(mut, c("WT", "WM", "MT")),
         mut = fct_recode(mut,
                          `WT+mutagen` = "WM",
                          mutator = "MT"))

rm(all_data)



# A function that generate integer breaks for ggplot
integer_breaks <- function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))

GraphGeneEA_with_p_count <- function(df, bg, ABX, gene, include.stop = FALSE) {
  if(include.stop == FALSE) {
    df <- filter(df, ACTION < 100)
    bg <- filter(bg, ACTION < 100)
  }
  plt_data <- df %>%
    filter(abx == ABX, GN == gene)
  bg_data <- bg %>%
    mutate(EA.bin = cut(ACTION, breaks = 10*0:10)) %>%
    mutate(EA.bin = as.numeric(EA.bin) * 10 -5) %>%
    filter(AF > 0.1) %>%
    select(-mut)
  bg_count <- bg_data %>%
    group_by(EA.bin) %>%
    summarise(count = n())
  mutColor <- data.frame(mut = c("WT", "WT+mutagen", "mutator"), color = c("black", "blue", "red"), stringsAsFactors = FALSE) %>%
    filter(mut %in% plt_data$mut)
  KS.p <- plt_data %>%
    group_by(mut) %>%
    summarise(ps = suppressWarnings(ks.test(ACTION, bg_data$ACTION,alternative = "less"))$p.value) %>%
    ungroup() %>%
    mutate(ps = signif(ps, digits = 2)) %>%
    mutate(mut_title = paste0(mut, "\n(EA-KS p = ", ps, ")")) %>%
    mutate(mut2 = factor(mut, labels = mut_title))
  mut_adj <- KS.p %>%
    select(mut, mut2)
  plt_data <- left_join(plt_data, KS.p, by = "mut") %>%
    mutate(mut = mut2)
  bg_factor <- nrow(bg)
  evolve_factor <- plt_data %>%
    group_by(mut) %>% 
    summarise(tot_count =  n()) %>%
    ungroup() %>%
    select(mut, tot_count) %>%
    mutate(scaling_factor = tot_count/bg_factor)
  bg_count <- crossing(bg_count, evolve_factor) %>%
    mutate(y = count*scaling_factor) %>%
    select(EA.bin, mut, y)
  plt <- ggplot() +
    geom_histogram(data = plt_data, aes(x = ACTION, y=..count.., fill=mut), color = "black",position="identity", alpha=1, breaks = 10*0:10) +
    geom_bar(data = bg_count, aes(x = EA.bin, y = y, group = mut), width = 10, alpha = 0.2, stat="identity") +
    theme(text = element_text(size=15)) +
    xlim(0,100) +
    scale_fill_manual(values = mutColor$color)+
    scale_y_continuous("Count",
                       breaks = integer_breaks) +
    labs(x = "EA bins") +
    theme_classic() +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 14), # set font for label bars
          axis.text = element_text(size = 12), # set font for axis numbers
          axis.title = element_text(size = 14), # set font for axis titles
          plot.title = element_text(size = 16, face="italic")) + 
    facet_wrap(~mut, scales = "free_y", nrow =1) +
    ggtitle(gene)
  return(plt)
}

GraphGeneSIFT_with_p_count <- function(df, bg, ABX, gene, include.stop = FALSE) {
  if(include.stop == FALSE) {
    df <- filter(df, ACTION < 100)
    bg <- filter(bg, ACTION < 100)
  }
  plt_data <- df %>%
    filter(abx == ABX, GN == gene)
  bg_data <- bg %>%
    mutate(SIFT.bin = cut(SIFT_adj, breaks = 10*0:10)) %>%
    mutate(SIFT.bin = as.numeric(SIFT.bin) * 10 -5) %>%
    filter(AF > 0.1) %>%
    select(-mut)
  bg_count <- bg_data %>%
    group_by(SIFT.bin) %>%
    summarise(count = n())
  mutColor <- data.frame(mut = c("WT", "WT+mutagen", "mutator"), color = c("black", "blue", "red"), stringsAsFactors = FALSE) %>%
    filter(mut %in% plt_data$mut)
  KS.p <- plt_data %>%
    group_by(mut) %>%
    summarise(ps = suppressWarnings(ks.test(SIFT_adj, bg_data$SIFT_adj,alternative = "less"))$p.value) %>%
    ungroup() %>%
    mutate(ps = signif(ps, digits = 2)) %>%
    mutate(mut_title = paste0(mut, "\n(SIFT_adj_KS p = ", ps, ")")) %>%
    mutate(mut2 = factor(mut, labels = mut_title))
  mut_adj <- KS.p %>%
    select(mut, mut2)
  plt_data <- left_join(plt_data, KS.p, by = "mut") %>%
    mutate(mut = mut2)
  bg_factor <- nrow(bg)
  evolve_factor <- plt_data %>%
    group_by(mut) %>% 
    summarise(tot_count =  n()) %>%
    ungroup() %>%
    select(mut, tot_count) %>%
    mutate(scaling_factor = tot_count/bg_factor)
  bg_count <- crossing(bg_count, evolve_factor) %>%
    mutate(y = count*scaling_factor) %>%
    select(SIFT.bin, mut, y)
  plt <- ggplot() +
    geom_histogram(data = plt_data, aes(x = SIFT_adj, y=..count.., fill=mut), color = "black",position="identity", alpha=1, breaks = 10*0:10) +
    geom_bar(data = bg_count, aes(x = SIFT.bin, y = y, group = mut), width = 10, alpha = 0.2, stat="identity") +
    theme(text = element_text(size=15)) +
    xlim(0,100) +
    scale_fill_manual(values = mutColor$color)+
    scale_y_continuous("Count",
                       breaks = integer_breaks) +
    labs(x = "SIFT_adj bins") +
    theme_classic() +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 14), # set font for label bars
          axis.text = element_text(size = 12), # set font for axis numbers
          axis.title = element_text(size = 14), # set font for axis titles
          plot.title = element_text(size = 16, face="italic")) + 
    facet_wrap(~mut, scales = "free_y", nrow =1) +
    ggtitle(gene)
  return(plt)
}


GraphGeneEA_with_p_count(evolve, random, ABX = "CIP", gene = "gyrA", include.stop = TRUE) 
GraphGeneSIFT_with_p_count(evolve, random, ABX = "CIP", gene = "gyrA", include.stop = TRUE) 

GraphGeneEA_with_p_count(evolve, random, ABX = "COL", gene = "basS", include.stop = TRUE) 
GraphGeneSIFT_with_p_count(evolve, random, ABX = "COL", gene = "basS", include.stop = TRUE) 

GraphGeneEA_with_p_count(evolve, random, ABX = "COL", gene = "ynjC", include.stop = TRUE) 
GraphGeneSIFT_with_p_count(evolve, random, ABX = "COL", gene = "ynjC", include.stop = TRUE) 
ggsave("plot/ALE_single_gene_EAdist/ynjC_SIFT.pdf", height = 3, width = 5, units = "in")

# ggplot(random) +
#   geom_histogram(aes(x = SIFT_adj, y=..count..), color = "black",position="identity", alpha=1, breaks = 10*0:10) +
#   ggtitle("Distribution of adjusted SIFT scores for random mutations")
#   

figures <- tibble(abx = c("CIP", "CIP", "COL", "COL"),
                  gene = c("gyrA", "parC","basS", "basR"),
                  width = c(6, 6, 6, 6)*1.2) %>%
  mutate(figs = map2(abx, gene, ~GraphGeneEA_with_p_count(evolve, random, ABX = .x, gene = .y, include.stop = TRUE)))
pmap(list(figures$figs, figures$gene, figures$width), ~ggsave(filename = paste0("plot/ALE_single_gene_EAdist/", ..2, ".pdf"), 
                                                              plot = ..1, height = 3, width = ..3, units = "in"))


ynjC <- evolve %>%
  filter(abx == "COL", GN == "ynjC")
