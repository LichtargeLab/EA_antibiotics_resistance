# Author: Chen Wang, Dec 2021
# Perform downsample analysis with the ALE dataset.
# The downsample analysis took days to run on a computing cluster.
# The results are included in the output directory.

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(stringr)
library(furrr)
library(scales)
library(ggh4x)
# plan(multicore, workers = 20)

############################################################################
# Functions
############################################################################
# KS and POI analyses
KS_Freq <- function(evolve, bg) {
  # Prepare gene length and gene name info
  len.map <- read_csv("data/len_map.csv",
                      col_types = "ccd")
  totallen <- sum(len.map$Len)
  
  avg.EA <- mean(bg$ACTION)
  output <- evolve %>%
    group_by(b_num) %>%
    summarize(ps = ks.test(ACTION, bg$ACTION,alternative = "less")$p.value,
              times = n(),
              strain_count = n_distinct(strain),
              sumEA = sum(ACTION)
    ) %>%
    left_join(len.map, by = "b_num") %>%
    mutate(expect_mut_count = nrow(evolve)/totallen*Len,
           p_poisson = ppois(times-1, expect_mut_count, lower.tail = FALSE),
           expect_sumEA = expect_mut_count * avg.EA,
           adj.sumEA = sumEA - expect_sumEA) %>%
    mutate(EAKS.rank = rank(ps),
           EAsum.rank = rank(desc(adj.sumEA)),
           Freq.rank = rank(p_poisson))
  return(output)
}

# Downsample function
DS <- function(evolve.df, bg.df, seed = 100) {
  unique.strains <- unique(evolve.df$strain)
  sampleN <- length(unique.strains)
  set.seed(seed)
  downSamp <- data.frame(sample_count = 1:(sampleN-1)) %>%
    crossing(rep = 1:10) %>% # for each sample size, try 10 times
    mutate(picks = map(sample_count, ~sample(unique.strains, .))) 
  downSamp.full <- tibble(sample_count = sampleN) %>%
    mutate(rep = 1,
           picks = list(unique.strains))
  downSamp <- bind_rows(downSamp, downSamp.full) %>%
    mutate(results = future_map(picks, ~KS_Freq(filter(evolve.df, strain %in% .), bg.df), .progress = TRUE))
  return(downSamp)
}


# load data
all_data <- read_csv("data/ALE_mutations_with_SIFT.csv",
                     col_types = "dccccddcccd")

random <- filter(all_data, mut == "random")


# founder mutations
d0 <- filter(all_data, abx == "d0")

evolve <- filter(all_data, abx %in% c("CIP", "COL")) %>%
  filter(AF > 0.1) %>%
  anti_join(., d0, by = c("POS", "REF", "ALT", "mut")) %>%
  group_by(mut, abx) %>%
  nest()

cip.MT.DS <- DS(evolve$data[[1]], bg.df = random, seed = 343)
cip.WM.DS <- DS(evolve$data[[3]], bg.df = random, seed = 23)
cip.WT.DS <- DS(evolve$data[[5]], bg.df = random, seed = 10)
col.MT.DS <- DS(evolve$data[[2]], bg.df = random, seed = 949)
col.WM.DS <- DS(evolve$data[[4]], bg.df = random, seed = 829)
col.WT.DS <- DS(evolve$data[[6]], bg.df = random, seed = 154)

# save(list = c("cip.MT.DS", "cip.WM.DS", "cip.WT.DS",
#               "col.MT.DS", "col.WM.DS", "col.WT.DS"),
#      file = "output/downsampling_results.RData")

load("output/downsampling_results.RData")

# Positive genes
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

# Plot the rankings of True positive genes
ExtractTP <- function(df, TP) {
  output <- df %>%
    filter(GENE_name %in% TP) %>%
    select(GENE_name, EAKS.rank, EAsum.rank, Freq.rank)
  return(output)
}


PlotDSTP <- function(df, TP, type = c("mean", "median"), scale_rank = FALSE) {
  if(scale_rank == TRUE) {
    workdf <- df %>%
      mutate(results = map(results, ~mutate(., across(contains(".rank"), function(x){x/max(x)}))))
  } else {
    workdf <- df
  }
  workdf <- workdf %>%
    mutate(data = map(results, ~ExtractTP(., TP = TP))) %>%
    select(sample_count, rep, data) %>%
    unnest(cols = c(data)) %>%
    mutate(GENE_name = factor(GENE_name, levels = TP)) %>%
    pivot_longer(cols = contains(".rank"), names_to = "method", values_to = "rank") %>%
    mutate(method = str_sub(method, 1, -6)) %>%
    group_by(sample_count, GENE_name, method) %>%
    summarize(rank.avg = mean(rank, na.rm = TRUE),
              rank.std = sd(rank, na.rm = TRUE),
              rank.median = median(rank, na.rm = TRUE),
              rank.qt75 = quantile(rank, 3/4, na.rm = TRUE),
              rank.qt25 = quantile(rank, 1/4, na.rm = TRUE),
              rank.min = min(rank, na.rm = TRUE),
              rank.max = max(rank, na.rm = TRUE)) %>%
    ungroup() 
  if (type == "mean") {
    output <- workdf %>%
      ggplot(aes(x = sample_count, y = rank.avg, color = method)) +
      geom_line() +
      geom_point(aes(shape = method), size = 2) +
      scale_y_reverse(breaks= pretty_breaks()) +
      xlab("sample size") +
      ylab("Mean rank") +
      facet_wrap(GENE_name ~., scales = "free_y", ncol = 2) +
      theme_bw() +
      theme(strip.text = element_text(size = 11)) +
      force_panelsizes(rows = unit(1.2, "in"), 
                       cols = unit(2.2, "in"), 
                       TRUE)
  } else if (type == "median") {
    output <- workdf %>%
      ggplot(aes(x = sample_count, y = rank.median, color = method)) +
      geom_line() +
      geom_point(aes(shape = method), size = 2) +
      geom_errorbar(aes(ymin = rank.qt25, ymax = rank.qt75),
                    width = 0.2, position = position_dodge(0.05)) +
      scale_y_reverse(breaks= pretty_breaks()) +
      xlab("sample size") +
      ylab("Median rank") +
      facet_wrap(GENE_name ~., scales = "free_y", ncol = 2) +
      theme_bw() +
      theme(strip.text = element_text(size = 11)) +
      force_panelsizes(rows = unit(1.2, "in"), 
                       cols = unit(2.2, "in"), 
                       TRUE)
  }
  
  return(output)
}


PlotDSTP(cip.WT.DS, cip.genes$GENE_name, type = "mean") +
  ggtitle("CIP - WT")
PlotDSTP(cip.WT.DS, cip.genes$GENE_name, type = "median") +
  ggtitle("CIP - WT") +
  theme(
        legend.text = element_text(size =12),
        strip.text.x = element_text(size = 14), # set font for label bars
        axis.text = element_text(size = 12), # set font for axis numbers
        axis.title = element_text(size = 14), # set font for axis titles
        title = element_text(size = 14))
ggsave("plot/downsampling/ALE_CIP_WT.jpeg", width = 7, height = 7, units = "in")

PlotDSTP(cip.WM.DS, cip.genes$GENE_name, type = "mean") +
  ggtitle("CIP - WT+mutagen")
PlotDSTP(cip.WM.DS, cip.genes$GENE_name, type = "median") +
  ggtitle("CIP - WT+mutagen") +
  theme(
    legend.text = element_text(size =12),
    strip.text.x = element_text(size = 14), # set font for label bars
    axis.text = element_text(size = 12), # set font for axis numbers
    axis.title = element_text(size = 14), # set font for axis titles
    title = element_text(size = 14))
ggsave("plot/downsampling/ALE_CIP_WM.jpeg", width = 7, height = 8, units = "in")


PlotDSTP(cip.MT.DS, cip.genes$GENE_name, type = "mean") +
  ggtitle("CIP - mutator")
PlotDSTP(cip.MT.DS, cip.genes$GENE_name, type = "median") +
  ggtitle("CIP - mutator") +
  theme(
    legend.text = element_text(size =12),
    strip.text.x = element_text(size = 14), # set font for label bars
    axis.text = element_text(size = 12), # set font for axis numbers
    axis.title = element_text(size = 14), # set font for axis titles
    title = element_text(size = 14))
ggsave("plot/downsampling/ALE_CIP_MT.jpeg", width = 7, height = 9, units = "in")


PlotDSTP(col.WT.DS, col.genes$GENE_name, type = "mean") +
  ggtitle("COL - WT")
PlotDSTP(col.WT.DS, col.genes$GENE_name, type = "median") +
  ggtitle("COL - WT") +
  theme(
    legend.text = element_text(size =12),
    strip.text.x = element_text(size = 14), # set font for label bars
    axis.text = element_text(size = 12), # set font for axis numbers
    axis.title = element_text(size = 14), # set font for axis titles
    title = element_text(size = 14))
ggsave("plot/downsampling/ALE_COL_WT.jpeg", width = 7, height = 7, units = "in")


PlotDSTP(col.WM.DS, col.genes$GENE_name, type = "mean") +
  ggtitle("COL - WT+mutagen")
PlotDSTP(col.WM.DS, col.genes$GENE_name, type = "median") +
  ggtitle("COL - WT+mutagen") +
  theme(
    legend.text = element_text(size =12),
    strip.text.x = element_text(size = 14), # set font for label bars
    axis.text = element_text(size = 12), # set font for axis numbers
    axis.title = element_text(size = 14), # set font for axis titles
    title = element_text(size = 14))
ggsave("plot/downsampling/ALE_COL_WM.jpeg", width = 7, height = 9, units = "in")


PlotDSTP(col.MT.DS, col.genes$GENE_name, type = "mean") +
  ggtitle("COL - mutator")
PlotDSTP(col.MT.DS, col.genes$GENE_name, type = "median") +
  ggtitle("COL - mutator") +
  theme(
    legend.text = element_text(size =12),
    strip.text.x = element_text(size = 14), # set font for label bars
    axis.text = element_text(size = 12), # set font for axis numbers
    axis.title = element_text(size = 14), # set font for axis titles
    title = element_text(size = 14))
ggsave("plot/downsampling/ALE_COL_MT.jpeg", width = 7, height = 9, units = "in")
