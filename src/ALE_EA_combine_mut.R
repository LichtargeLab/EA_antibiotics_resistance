library(tidyverse)


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


############################################################################
# Prepare inputs and bg distribution
############################################################################

all_data <- read_csv("data/ALE_mutations_with_SIFT.csv",
                     col_types = "dccccddcccd")

# randomly simulated mutatons
random <- filter(all_data, mut == "random")

# founder mutations
d0 <- filter(all_data, abx %in% c("d0", "passage"))

# Founder mutations include all the mutations that are reported by breseq (AF 0.05)
# But in the testing, only mutations with AF > 0.1 are used.
evolve <- filter(all_data, abx %in% c("CIP", "COL")) %>%
  filter(AF > 0.1) %>%
  anti_join(., d0, by = c("POS", "REF", "ALT")) %>%
  group_by(abx) %>%
  nest()

evolve_analysis <- evolve %>%
  mutate(result = map(data, ~KS_Freq(., random)))

cip_result <- evolve_analysis$result[[1]]
col_result <- evolve_analysis$result[[2]]

openxlsx::write.xlsx(list(cip_result, col_result), 
                     sheetName = c("cipro", "colistin"),
                     "output/ALE_EA_combine_mut_output.xlsx")
