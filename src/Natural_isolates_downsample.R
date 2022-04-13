# Downsampling assay with known drivers.

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(ggplot2)
library(furrr)
library(scales)
library(cowplot)
library(plotROC)
plan(multicore, workers = 20)

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

ClinicalSumEA <- function(SNP_df) {
  count_df <- SNP_df %>%
    select(strain, resis) %>%
    unique %>%
    group_by(resis) %>%
    summarise(count = n())
  unique_mut <- SNP_df %>%
    select(b_num, GENE_name, SUB, ACTION, resis) %>%
    unique() %>%
    arrange(resis) %>%
    group_by(b_num, GENE_name, SUB, ACTION) %>%
    summarize(unique = paste0(resis, collapse = ""))
  sumEA_output <- SNP_df %>%
    left_join(unique_mut) %>%
    filter(unique %in% c("R", "S")) %>%
    group_by(b_num, GENE_name, resis) %>%
    summarise(sumEA = sum(ACTION)) %>%
    left_join(count_df) %>%
    mutate(avg.sumEA = sumEA/count) %>%
    select(-count, -sumEA) %>%
    pivot_wider(names_from = resis, values_from = avg.sumEA, values_fill = 0) %>%
    mutate(diff = R - S) %>%
    filter(R > 0) %>%
    arrange(desc(diff)) %>%
    select(b_num, GENE_name, R, S, diff) %>%
    rename(sumEA.per.R.strain = R,
           sumEA.per.S.strain = S) %>%
    ungroup() %>%
    mutate(sumEA.R.rank = rank(desc(sumEA.per.R.strain)),
           sumEA.diff.rank = rank(desc(diff))) 
  return(sumEA_output)
}


# Downsample function for natural isolates data
DS_NI <- function(mut.df, seed = 100) {
  mut.R <- mut.df %>%
    filter(resis == "R")
  mut.S <- mut.df %>%
    filter(resis == "S")
  unique.strains.R <- unique(mut.R$strain)
  sampleN <- length(unique.strains.R)
  set.seed(seed)
  downSamp <- data.frame(sample_count = seq(1, sampleN%/%10) * 10) %>%
    crossing(rep = 1:10) %>% # for each sample size, try 10 times
    mutate(picks = map(sample_count, ~sample(unique.strains.R, .))) 
  downSamp.full <- tibble(sample_count = sampleN) %>%
    mutate(rep = 1,
           picks = list(unique.strains.R))
  downSamp <- bind_rows(downSamp, downSamp.full)
  # KS and Freq tests
  mut.bg <- mut.S %>%
    select(-strain) %>%
    unique()
  mut.R.filt <- mut.R %>%
    anti_join(mut.S, by = c("b_num", "SUB"))
  downSamp <- downSamp %>%
    mutate(KSFreq_output = future_map(picks, ~KS_Freq(filter(mut.R.filt, strain %in% .), mut.bg), .progress = TRUE))
  # EA sum test
  downSamp <- downSamp %>%
    mutate(EAsum_output = future_map(picks, ~ClinicalSumEA(bind_rows(filter(mut.R, strain %in% .), mut.S)), .progress = TRUE))
  # Combine results
  CombineResults <- function(KSFreq.df, EAsum.df) {
    output <- select(KSFreq.df, b_num, GENE_name, Len, times, strain_count, EAKS.rank, Freq.rank) %>%
      left_join(select(EAsum.df, b_num, EAsum.rank = sumEA.diff.rank)) %>%
      relocate(b_num, GENE_name, Len, times, strain_count, EAKS.rank, EAsum.rank, Freq.rank)
    return(output)
  }
  downSamp <- downSamp %>%
    mutate(results = map2(KSFreq_output, EAsum_output, CombineResults)) %>%
    select(-KSFreq_output, -EAsum_output)
  return(downSamp)
}



cip <- readRDS("data/Natural_isolates_cipro_mutations_with_SIFT.RDS")

cip.NI.DS <- DS_NI(cip, seed = 230)

col <- readRDS("data/Natural_isolates_colistin_mutations_with_SIFT.RDS")

col.NI.DS <- DS_NI(col, seed = 304)

# save(list = c("cip.NI.DS", "col.NI.DS"),
#      file = "output/Natural_isolates_downsampling_results.RData")

load("output/Natural_isolates_downsampling_results.RData")


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


ExtractTP <- function(df, TP) {
  output <- df %>%
    filter(GENE_name %in% TP) %>%
    select(GENE_name, EAKS.rank, EAsum.rank, Freq.rank)
  return(output)
}

library(ggh4x)

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

PlotDSTP(cip.NI.DS, c("gyrA", "parC", "parE", "acrR"), type = "median") +
  theme(legend.position="bottom",
        legend.text = element_text(size =12),
        strip.text.x = element_text(size = 14, face = "italic"), # set font for label bars
        axis.text = element_text(size = 12), # set font for axis numbers
        axis.title = element_text(size = 14), # set font for axis titles
        title = element_text(size = 14))
ggsave("plot/downsampling/Natural_isolates_cipro.pdf", width = 7, height = 5, units = "in")


PlotDSTP(col.NI.DS, c("basS", "basR"), type = "median") +
  theme(legend.position="bottom",
        legend.text = element_text(size =12),
        strip.text.x = element_text(size = 14, face = "italic"), # set font for label bars
        axis.text = element_text(size = 12), # set font for axis numbers
        axis.title = element_text(size = 14), # set font for axis titles
        title = element_text(size = 14))
ggsave("plot/downsampling/Natural_isolates_colistin.pdf", width = 7, height = 3, units = "in")




cip.NI.DS %>%
  mutate(mutated_gene_count = map_dbl(results, nrow)) %>%
  select(sample_count, rep, mutated_gene_count) %>%
  mutate(sample_count = as.factor(sample_count)) %>%
  ggplot(aes(x = sample_count, y = mutated_gene_count)) + #as.factor function is important here, otherwise only one box for each fill
  geom_boxplot(notch=FALSE, outlier.shape=NA, alpha=0.5) +
  xlab("sample_count") +
  scale_y_reverse(breaks= pretty_breaks()) +
  ggtitle("CIP-clinical-all genes") +
  theme_bw()

col.NI.DS %>%
  mutate(mutated_gene_count = map_dbl(results, nrow)) %>%
  select(sample_count, rep, mutated_gene_count) %>%
  mutate(sample_count = as.factor(sample_count)) %>%
  ggplot(aes(x = sample_count, y = mutated_gene_count)) + #as.factor function is important here, otherwise only one box for each fill
  geom_boxplot(notch=FALSE, outlier.shape=NA, alpha=0.5) +
  xlab("sample_count") +
  scale_y_reverse(breaks= pretty_breaks()) +
  ggtitle("COL-clinical-all genes") +
  theme_bw()


cipro_TP_output <- cip.NI.DS %>%
  mutate(data = map(results, ~ExtractTP(., TP = c("gyrA", "parC", "parE", "acrR")))) %>%
  select(sample_count, rep, data) %>%
  unnest(cols = c(data))

colistin_TP_output <- col.NI.DS %>%
  mutate(data = map(results, ~ExtractTP(., TP = c("basS", "basR")))) %>%
  select(sample_count, rep, data) %>%
  unnest(cols = c(data))

openxlsx::write.xlsx(list(cipro_TP_output, colistin_TP_output),
                     "output/Natural_isolates_downsampling_TP.xlsx",
                     sheetName = c("cipro", "colistin"))
