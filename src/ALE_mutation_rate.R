library(tidyverse)
library(magrittr)
library(ggpubr)
library(rstatix)

all_data <- read_csv("data/ALE_mutations_with_SIFT.csv")

d0 <- filter(all_data, abx == "d0")

mutation_rate <- all_data %>%
  filter(abx %in% c("CIP", "COL")) %>%
  filter(AF > 0.1) %>%
  anti_join(., d0, by = c("POS", "REF", "ALT", "mut")) %>%
  group_by(abx, mut, strain) %>%
  summarize(coding_mut_count = n()) %>%
  ungroup()

days <- select(mutation_rate, abx, mut) %>%
  unique() %>%
  mutate(day = c(14, 16, 16, 21,17,17))

mutation_rate <- mutation_rate %>%
  left_join(days) %>%
  mutate(coding_mut_per_day = coding_mut_count/day) %>%
  mutate(mut = factor(mut, levels = c("WT", "WM", "MT")),
         mut = fct_recode(mut, WT = "WT",
                          `WT+mutagen` = "WM",
                          mutator = "MT"),
         abx = fct_recode(abx, ciprofloxacin = "CIP",
                          colistin = "COL"))

write_csv(mutation_rate, "output/ALE_mutation_rate_per_sample.csv")

ggplot(mutation_rate) +
  geom_dotplot(aes(x = mut, y = coding_mut_per_day, fill = mut), binaxis = "y", stackdir = "center",
               position=position_dodge(1), binwidth = 0.04, dotsize = 1.5) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  scale_fill_manual(values = c("black", "blue", "red")) +
  theme_classic() +
  theme(legend.position = "top",
        legend.justification='left',
        axis.text = element_text(size = 18, face = "bold"), # set font for axis numbers
        legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.8),
        axis.ticks = element_line(colour = "black", size = 0.8),
        axis.title = element_text(size = 18, face = "bold")) + # set font for axis titles
  xlab(NULL) +
  ylab("Average coding mutation count per day") +
  facet_grid(.~abx)


my_comparisons <- list(c("WT", "WT+mutagen"), c("WT+mutagen", "mutator"), c("WT", "mutator"))

add_logticks  <- function (base = 10, sides = "bl", scaled = TRUE, 
                           short = unit(0.1, "cm"), mid = unit(0.2, "cm"),  long = unit(0.3, "cm"), 
                           colour = "black",  size = 0.5, linetype = 1, alpha = 1, color = NULL, 
                           data =data.frame(x = NA),... )   {
  if (!is.null(color)) 
    colour <- color
  layer(geom = "logticks", params = list(base = base, 
                                         sides = sides, scaled = scaled, short = short, 
                                         mid = mid, long = long, colour = colour, size = size, 
                                         linetype = linetype, alpha = alpha, ...), 
        stat = "identity", data =data , mapping = NULL, inherit.aes = FALSE, 
        position = "identity",
        show.legend = FALSE)
}


pairwise.test <- mutation_rate %>%
  group_by(abx) %>%
  t_test(coding_mut_per_day ~ mut, p.adjust.method = 'bonferroni') %>%
  adjust_pvalue(method = 'bonferroni') %>%
  add_significance() %>%
  add_xy_position(x = "mut") %>%
  mutate(y.position = log10(c(100, 500, 200, 100, 500, 200))) %>%
  mutate(p.label = paste0("p=", signif(p.adj, 2)))


mutation_rate %>%
  ggdotplot(x = "mut", y = "coding_mut_per_day",
            fill = "mut", palette = c("black", "blue", "red"),
            facet.by = "abx", binwidth = 0.04, dotsize = 2) +
  scale_y_log10() +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=8.5),
        axis.text.x = element_text(angle = 90, size=8.5), # set font for axis numbers
        axis.text.y = element_text(size=8.5), 
        axis.title = element_text(size=8.5),
        strip.text = element_text(size=8.5),
        legend.title = element_blank()
  ) + # set font for axis titles
  xlab(NULL) +
  ylab("Average coding mutation count per day") +
  stat_pvalue_manual(
    pairwise.test,  label = "p.label", tip.length = 0, size = 2.4
  ) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  add_logticks(sides = "l", data = data.frame(x= NA, abx = "ciprofloxacin"))

ggsave("plot/ALE_mutation_rate.pdf", width = 2.4, height = 3.3, units = "in")


