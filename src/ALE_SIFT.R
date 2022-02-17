# Author: Chen Wang, Dec 2021
# Replace EA scores with adjusted SIFT scores in EAKS and EAsum, and repeat the analyses 
# on ALE dataset

library(tidyverse)
library(RobustRankAggreg)

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

# Prepare gene list for aggregation. Remove NAs.
NameRankings <- function(df, var) {
  workdf <- select(df, b_num, all_of(var)) %>%
    filter(!is.na(!!sym(var))) %>%
    arrange(!!sym(var))
  return(workdf$b_num)
}
# Combines gene rankings from three mutation loads in the output dataframe
CombineRankings_SIFT <- function(df) {
  SIFT_adj_KS_combine <- aggregateRanks(list(NameRankings(df, "SIFT_adj_KS.rank_MT"),
                                             NameRankings(df, "SIFT_adj_KS.rank_WM"),
                                             NameRankings(df, "SIFT_adj_KS.rank_WT"))) %>%
    mutate(SIFT_adj_KS_combine = rank(Score)) %>%
    select(SIFT_adj_KS_score = Score, b_num = Name, SIFT_adj_KS_combine)
  SIFT_adj_sum_combine <- aggregateRanks(list(NameRankings(df, "SIFT_adj_sum.rank_MT"),
                                              NameRankings(df, "SIFT_adj_sum.rank_WM"),
                                              NameRankings(df, "SIFT_adj_sum.rank_WT"))) %>%
    mutate(SIFT_adj_sum_combine = rank(Score)) %>%
    select(SIFT_adj_sum_score = Score, b_num = Name, SIFT_adj_sum_combine)
  Freq_combine <- aggregateRanks(list(NameRankings(df, "Freq.rank_MT"),
                                      NameRankings(df, "Freq.rank_WM"),
                                      NameRankings(df, "Freq.rank_WT"))) %>%
    mutate(Freq_combine = rank(Score)) %>%
    select(Freq_score = Score, b_num = Name, Freq_combine)
  output <- df %>%
    left_join(SIFT_adj_KS_combine, by = "b_num") %>%
    left_join(SIFT_adj_sum_combine, by = "b_num") %>%
    left_join(Freq_combine, by = "b_num") %>%
    arrange(SIFT_adj_KS_combine)
  return(output)
}

############################################################################
# Prepare inputs and bg distribution
############################################################################

# In the dataset, SIFT scores are adjusted with 100*(1-SIFT), so that it has
# the same scale as EA (ACTION). With 100 being most deleterious, and 0 being
# the lest impactful mutations
all_data <- read_csv("data/ALE_mutations_with_SIFT.csv",
                     col_types = "dccccddcccd")

# Some mutations doesn't have SIFT scores annotated, so they are removed from
# the analysis.
all_data <- all_data %>%
  filter(!is.na(SIFT_adj))

# randomly simulated mutations
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
  mutate(result = map(data, ~KS_Freq_SIFT(., random)))

saveRDS(evolve_analysis, "output/ALE_SIFT_details.RDS")

evolve_output <- evolve_analysis %>%
  unnest(cols = c(result)) %>%
  ungroup() %>%
  select(mut, abx, b_num, GENE_name, SIFT_adj_KS.rank, SIFT_adj_sum.rank, Freq.rank) 
cipro <- evolve_output %>%
  filter(abx == "CIP") %>%
  select(-abx) %>%
  pivot_wider(id_cols = c(b_num, GENE_name), names_from = c(mut), values_from = c(SIFT_adj_KS.rank, SIFT_adj_sum.rank, Freq.rank)) %>%
  CombineRankings_SIFT()
colistin <- evolve_output %>%
  filter(abx == "COL") %>%
  select(-abx) %>%
  pivot_wider(id_cols = c(b_num, GENE_name), names_from = c(mut), values_from = c(SIFT_adj_KS.rank, SIFT_adj_sum.rank, Freq.rank)) %>%
  CombineRankings_SIFT()



output_list <- list(cipro, colistin)


openxlsx::write.xlsx(list(cipro, colistin), 
                     sheetName = c("cipro", "colistin"),
                     "output/ALE_SIFT_output.xlsx")
