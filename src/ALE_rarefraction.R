library(tidyverse)
library(scales)
library(cowplot)

GetUniqueMut <- function(evolve.df, strains, genes = NULL) {
  if (!is.null(genes)) {
    evolve.df <- evolve.df %>%
      filter(GENE_name %in% genes)
  }
  mut.df <- evolve.df %>%
    filter(strain %in% strains) %>%
    select(POS, REF, ALT, b_num, SUB) %>%
    unique()
  gene.df <- mut.df %>%
    select(b_num) %>%
    unique()
  output <- tibble(mutation_count = nrow(mut.df),
                   mutated_gene_count = nrow(gene.df))
  return(output)
}

RareFraction <- function(evolve.df, seed = 100, genes = NULL) {
  unique.strains <- unique(evolve.df$strain)
  sampleN <- length(unique.strains)
  set.seed(seed)
  rarefrac <- data.frame(sample_count = 1:(sampleN-1)) %>%
    mutate(rep = list(1:10)) %>%  # for each sample size, try 10 times
    unnest(cols = c(rep)) %>%
    mutate(picks = map(sample_count, ~sample(unique.strains, .))) 
  rarefrac.full <- tibble(sample_count = sampleN) %>%
    mutate(rep = 1,
           picks = I(list(unique.strains)))
  rarefrac <- bind_rows(rarefrac, rarefrac.full) %>%
    mutate(results = map(picks, ~GetUniqueMut(evolve.df, ., genes = genes))) %>%
    unnest(cols = c(results))
  return(rarefrac)
}

PlotRareFraction <- function(df, y_var, title) {
  output <- df %>% 
    mutate(sample_count = as.factor(sample_count)) %>%
    ggplot(aes_string(x = "sample_count", y = y_var)) + #as.factor function is important here, otherwise only one box for each fill
    geom_boxplot(notch=FALSE, outlier.shape=NA, alpha=0.5) +
    xlab("sample count") +
    scale_y_continuous(breaks= pretty_breaks()) +
    # ggtitle(title) +
    theme_bw() +
    theme(
      legend.text = element_text(size =12),
      strip.text.x = element_text(size = 14), # set font for label bars
      axis.text = element_text(size = 10), # set font for axis numbers
      axis.title = element_text(size = 14), # set font for axis titles
      title = element_text(size = 14))
  return(output)
}


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

# load data
all_data <- read_csv("data/ALE_mutations_with_SIFT.csv",
                     col_types = "dccccddcccd")
len.map <- read_csv("data/len_map.csv")

# founder mutations
d0 <- filter(all_data, abx == "d0")

evolve <- filter(all_data, abx %in% c("CIP", "COL")) %>%
  filter(AF > 0.1) %>%
  anti_join(., d0, by = c("POS", "REF", "ALT", "mut")) %>%
  left_join(len.map) %>%
  group_by(mut, abx) %>%
  nest() %>%
  mutate(TP = ifelse(abx == "CIP",
                     list(cip.genes$GENE_name),
                     list(col.genes$GENE_name)))

rarefraction.output <- evolve %>%
  ungroup %>%
  mutate(seed = c(343,949,23,829,10,154)) %>%
  mutate(result_all = map2(data, seed, ~RareFraction(evolve.df = .x, seed = .y)),
         result_TP = pmap(list(data, seed, TP), ~RareFraction(evolve.df = ..1, seed = ..2, genes = ..3))) %>%
  mutate(
    #all_gene_plot = pmap(list(result_all, abx, mut), ~PlotRareFraction(..1, "mutated_gene_count", paste0(..2, "-", ..3, "-all genes"))),
    #all_mut_plot = pmap(list(result_all, abx, mut), ~PlotRareFraction(..1, "mutation_count", paste0(..2, "-", ..3, "-all genes"))),
    TP_gene_plot = pmap(list(result_TP, abx, mut), ~PlotRareFraction(..1, "mutated_gene_count", paste0(..2, "-", ..3, "-TP genes"))),
    # TP_mut_plot = pmap(list(result_TP, abx, mut), ~PlotRareFraction(..1, "mutation_count", paste0(..2, "-", ..3, "-TP genes")))
  ) %>%
  arrange(abx, desc(mut))

# plot_grid(plotlist = rarefraction.output$all_gene_plot, align = "hv", nrow = 2)
# plot_grid(plotlist = rarefraction.output$all_mut_plot, align = "hv", nrow = 2)
plot_grid(plotlist = rarefraction.output$TP_gene_plot, align = "hv", nrow = 2)
ggsave("plot/downsampling/ALE_TP.pdf", width = 12, height = 6, units = "in")
# plot_grid(plotlist = rarefraction.output$TP_mut_plot, align = "hv", nrow = 2)

output <- rarefraction.output %>%
  select(mut, abx, result_TP) %>%
  unnest(cols = c(result_TP)) %>%
  select(mut, abx, sample_count, rep, mutated_TP_gene = mutated_gene_count)

write_csv(output, "output/ALE_rarefraction.csv")
