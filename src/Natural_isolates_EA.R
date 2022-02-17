# Author: Chen Wang, Dec 2021
# Run natural isolate dataset with EA and frequency analyses
library(tidyverse)

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

############################################################################
# cipro
############################################################################
cip <- readRDS("data/Natrual_isolates_cipro_mutations_with_SIFT.RDS")

cip_KS_bg <- cip %>%
  filter(resis == "S") %>%
  select(-strain) %>%
  unique()

cip_KS_R <- cip %>%
  filter(resis == "R") %>%
  anti_join(cip_KS_bg, by = c("b_num", "SUB"))


# Note that the clinical EAsum has a slightly different approach comparing to the
# one implemented in KS_Freq function.
cip_KS_Freq <- KS_Freq(cip_KS_R, cip_KS_bg)
cip_EAsum <- ClinicalSumEA(cip)

saveRDS(cip_KS_Freq, "output/Natrual_isolates_cipro_EAKS_details.RDS")
saveRDS(cip_EAsum, "output/Natrual_isolates_cipro_EAsum_details.RDS")

cip_output <- select(cip_KS_Freq, b_num, GENE_name, EAKS.rank, Freq.rank) %>%
  full_join(select(cip_EAsum, b_num, EAsum.rank = sumEA.diff.rank)) %>%
  arrange(EAKS.rank)

rm(list = c("cip", "cip_KS_bg", "cip_KS_R"))

############################################################################
# colistin
############################################################################
col <- readRDS("data/Natrual_isolates_colistin_mutations_with_SIFT.RDS")

col_KS_bg <- col %>%
  filter(resis == "S") %>%
  select(-strain) %>%
  unique()

col_KS_R <- col %>%
  filter(resis == "R") %>%
  anti_join(col_KS_bg, by = c("b_num", "SUB"))


# Note that the clinical EAsum has a slightly different approach comparing to the
# one implemented in KS_Freq function.
col_KS_Freq <- KS_Freq(col_KS_R, col_KS_bg)
col_EAsum <- ClinicalSumEA(col)

saveRDS(col_KS_Freq, "output/Natrual_isolates_colistin_EAKS_details.RDS")
saveRDS(col_EAsum, "output/Natrual_isolates_colistin_EAsum_details.RDS")

col_output <- select(col_KS_Freq, b_num, GENE_name, EAKS.rank, Freq.rank) %>%
  full_join(select(col_EAsum, b_num, EAsum.rank = sumEA.diff.rank)) %>%
  arrange(EAKS.rank)


openxlsx::write.xlsx(list(cip_output, col_output), "output/Natrual_isolates_EA_output.xlsx",
                     sheetName = c("cipro", "colistin"))
