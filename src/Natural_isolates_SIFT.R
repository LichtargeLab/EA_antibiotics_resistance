# Author: Chen Wang, Dec 2021
# Replace EA scores with adjusted SIFT scores in EAKS and EAsum, and repeat the analyses 
# on natural isolate dataset.
library(tidyverse)

############################################################################
# Functions
############################################################################
# KS and POI analyses
KS_Freq_SIFT <- function(evolve, bg) {
  # Prepare gene length and gene name info
  len.map <- read_csv("data/len_map.csv")
  totallen <- sum(len.map$Len)
  
  avg.SIFT_adj <- mean(bg$SIFT_adj)
  output <- evolve %>%
    group_by(b_num) %>%
    summarize(ps = ks.test(SIFT_adj, bg$SIFT_adj,alternative = "less")$p.value,
              times = n(),
              strain_count = n_distinct(strain),
              sumSIFT_adj = sum(SIFT_adj)
    ) %>%
    left_join(len.map, by = "b_num") %>%
    mutate(expect_mut_count = nrow(evolve)/totallen*Len,
           p_poisson = ppois(times-1, expect_mut_count, lower.tail = FALSE),
           expect_sumSIFT_adj = expect_mut_count * avg.SIFT_adj,
           adj.sumSIFT_adj = sumSIFT_adj - expect_sumSIFT_adj) %>%
    mutate(SIFT_adj_KS.rank = rank(ps),
           SIFT_adj_sum.rank = rank(desc(adj.sumSIFT_adj)),
           Freq.rank = rank(p_poisson))
  return(output)
}


ClinicalSumSIFT <- function(SNP_df) {
  count_df <- SNP_df %>%
    select(strain, resis) %>%
    unique %>%
    group_by(resis) %>%
    summarise(count = n())
  unique_mut <- SNP_df %>%
    select(b_num, GENE_name, SUB, SIFT_adj, resis) %>%
    unique() %>%
    arrange(resis) %>%
    group_by(b_num, GENE_name, SUB, SIFT_adj) %>%
    summarize(unique = paste0(resis, collapse = ""))
  sumSIFT_output <- SNP_df %>%
    left_join(unique_mut) %>%
    filter(unique %in% c("R", "S")) %>%
    group_by(b_num, GENE_name, resis) %>%
    summarise(sumSIFT = sum(SIFT_adj)) %>%
    left_join(count_df) %>%
    mutate(avg.sumSIFT = sumSIFT/count) %>%
    select(-count, -sumSIFT) %>%
    pivot_wider(names_from = resis, values_from = avg.sumSIFT, values_fill = 0) %>%
    mutate(diff = R - S) %>%
    filter(R > 0) %>%
    arrange(desc(diff)) %>%
    select(b_num, GENE_name, R, S, diff) %>%
    rename(sumSIFT.per.R.strain = R,
           sumSIFT.per.S.strain = S) %>%
    ungroup() %>%
    mutate(sumSIFT.R.rank = rank(desc(sumSIFT.per.R.strain)),
           sumSIFT.diff.rank = rank(desc(diff))) 
  return(sumSIFT_output)
}

############################################################################
# cipro
############################################################################
cip <- readRDS("data/Natural_isolates_cipro_mutations_with_SIFT.RDS") %>%
  filter(!is.na(SIFT_adj))

cip_KS_bg <- cip %>%
  filter(resis == "S") %>%
  select(-strain) %>%
  unique()

cip_KS_R <- cip %>%
  filter(resis == "R") %>%
  anti_join(cip_KS_bg, by = c("b_num", "SUB"))


# Note that the clinical EAsum has a slightly different approach comparing to the
# one implemented in KS_Freq function.
cip_KS_Freq <- KS_Freq_SIFT(cip_KS_R, cip_KS_bg)
cip_SIFTsum <- ClinicalSumSIFT(cip)

saveRDS(cip_KS_Freq, "output/Natural_isolates_cipro_SIFT_KS_details.RDS")
saveRDS(cip_SIFTsum, "output/Natural_isolates_cipro_SIFTsum_details.RDS")

cip_output <- select(cip_KS_Freq, b_num, GENE_name, SIFT_adj_KS.rank, Freq.rank) %>%
  full_join(select(cip_SIFTsum, b_num, SIFT_adj_sum.rank = sumSIFT.diff.rank)) %>%
  arrange(SIFT_adj_KS.rank)

rm(list = c("cip", "cip_KS_bg", "cip_KS_R"))

############################################################################
# colistin
############################################################################
col <- readRDS("data/Natural_isolates_colistin_mutations_with_SIFT.RDS") %>%
  filter(!is.na(SIFT_adj))

col_KS_bg <- col %>%
  filter(resis == "S") %>%
  select(-strain) %>%
  unique()

col_KS_R <- col %>%
  filter(resis == "R") %>%
  anti_join(col_KS_bg, by = c("b_num", "SUB"))


# Note that the clinical EAsum has a slightly different approach comparing to the
# one implemented in KS_Freq function.
col_KS_Freq <- KS_Freq_SIFT(col_KS_R, col_KS_bg)
col_SIFTsum <- ClinicalSumSIFT(col)

saveRDS(col_KS_Freq, "output/Natural_isolates_colistin_SIFT_KS_details.RDS")
saveRDS(col_SIFTsum, "output/Natural_isolates_colistin_SIFTsum_details.RDS")

col_output <- select(col_KS_Freq, b_num, GENE_name, SIFT_adj_KS.rank, Freq.rank) %>%
  full_join(select(col_SIFTsum, b_num, SIFT_adj_sum.rank = sumSIFT.diff.rank)) %>%
  arrange(SIFT_adj_KS.rank)


openxlsx::write.xlsx(list(cip_output, col_output), "output/Natural_isolates_SIFT_output.xlsx",
                     sheetName = c("cipro", "colistin"))
