# Author: Chen Wang, Dec 2021
# Run EAKS, EAsum and frequency based analyses on ALE dataset
library(tidyverse)
library(RobustRankAggreg)


############################################################################
# Functions
############################################################################
# KS and POI analyses
KS_Freq <- function(evolve, bg) {
  # Prepare gene length and gene name info
  len.map <- read_csv("data/len_map.csv")
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

# Prepare gene list for aggregation. Remove NAs.
NameRankings <- function(df, var) {
  workdf <- select(df, b_num, all_of(var)) %>%
    filter(!is.na(!!sym(var))) %>%
    arrange(!!sym(var))
  return(workdf$b_num)
}
# Combines gene rankings from three mutation loads in the output dataframe
CombineRankings <- function(df) {
  EAKS_combine <- aggregateRanks(list(NameRankings(df, "EAKS.rank_MT"),
                                    NameRankings(df, "EAKS.rank_WM"),
                                    NameRankings(df, "EAKS.rank_WT"))) %>%
    mutate(EAKS_combine = rank(Score)) %>%
    select(EAKS_score = Score, b_num = Name, EAKS_combine)
  EAsum_combine <- aggregateRanks(list(NameRankings(df, "EAsum.rank_MT"),
                                       NameRankings(df, "EAsum.rank_WM"),
                                       NameRankings(df, "EAsum.rank_WT"))) %>%
    mutate(EAsum_combine = rank(Score)) %>%
    select(EAsum_score = Score, b_num = Name, EAsum_combine)
  Freq_combine <- aggregateRanks(list(NameRankings(df, "Freq.rank_MT"),
                                      NameRankings(df, "Freq.rank_WM"),
                                      NameRankings(df, "Freq.rank_WT"))) %>%
    mutate(Freq_combine = rank(Score)) %>%
    select(Freq_score = Score, b_num = Name, Freq_combine)
  output <- df %>%
    left_join(EAKS_combine, by = "b_num") %>%
    left_join(EAsum_combine, by = "b_num") %>%
    left_join(Freq_combine, by = "b_num") %>%
    arrange(EAKS_combine)
  return(output)
}

############################################################################
# Prepare inputs and bg distribution
############################################################################

all_data <- read_csv("data/ALE_mutations_with_SIFT.csv",
                          col_types = "dccccddcccd")

# randomly simulated mutatons
random <- filter(all_data, mut == "random")

# founder mutations
d0 <- filter(all_data, abx == "d0")

# Founder mutations include all the mutations that are reported by breseq (AF 0.05)
# But in the testing, only mutations with AF > 0.1 are used.
evolve <- filter(all_data, abx %in% c("CIP", "COL")) %>%
  filter(AF > 0.1) %>%
  anti_join(., d0, by = c("POS", "REF", "ALT", "mut")) %>%
  group_by(mut, abx) %>%
  nest()




############################################################################
# KS and POI analyses with all mutations with AF above cutoff
############################################################################
evolve_analysis <- evolve %>%
  mutate(result = map(data, ~KS_Freq(., random)))

saveRDS(evolve_analysis, "output/ALE_EA_details.RDS")


evolve_output <- evolve_analysis %>%
  unnest(cols = c(result)) %>%
  ungroup() %>%
  select(mut, abx, b_num, GENE_name, EAKS.rank, EAsum.rank, Freq.rank) 
cipro <- evolve_output %>%
  filter(abx == "CIP") %>%
  select(-abx) %>%
  pivot_wider(id_cols = c(b_num, GENE_name), names_from = c(mut), values_from = c(EAKS.rank, EAsum.rank, Freq.rank)) %>%
  CombineRankings()
colistin <- evolve_output %>%
  filter(abx == "COL") %>%
  select(-abx) %>%
  pivot_wider(id_cols = c(b_num, GENE_name), names_from = c(mut), values_from = c(EAKS.rank, EAsum.rank, Freq.rank)) %>%
  CombineRankings()



output_list <- list(cipro, colistin)


openxlsx::write.xlsx(list(cipro, colistin), 
                     sheetName = c("cipro", "colistin"),
                     "output/ALE_EA_output.xlsx")


