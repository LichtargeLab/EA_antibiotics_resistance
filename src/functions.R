# Author: Chen Wang, Dec 2021

# Some functions


library(scales)

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}


GraphRankings <- function(graph.df, highlight.genes, xvar, yvar, title, 
                          xlab = "Frequency rank", ylab = "EA rank") {
  highlight.df <- graph.df %>%
    left_join(highlight.genes) %>%
    filter(!is.na(class))
  output.plot <- ggplot() +
    geom_point(data = filter(graph.df, !GENE_name %in% highlight.genes$GENE_name), aes(x = !!sym(xvar), y = !!sym(yvar)), fill = "gray60", color = "gray60", pch=21, alpha = 0.5) +
    geom_point(data = highlight.df, aes(x = !!sym(xvar), y = !!sym(yvar), fill = class), color = "black", pch=21, size = 2) +
    ggrepel::geom_text_repel(data = highlight.df, aes(x = !!sym(xvar), y = !!sym(yvar), label = GENE_name), nudge_x = -0.2, nudge_y = 0, fontface = "bold", size = 4) +
    geom_function(fun = function(x) x, color = "red", linetype = 2) +
    # scale_fill_manual(values = c("blue", "red")) +
    guides(fill=guide_legend(title="")) +
    scale_x_continuous(trans=reverselog_trans(10)) +
    scale_y_continuous(trans=reverselog_trans(10)) +
    labs(x = xlab, y = ylab) +
    coord_fixed() +
    theme_classic() +
    theme(
      title = element_text(face = "bold", size = 15),
      axis.title = element_text(face="bold", size = 15),
      axis.text = element_text(face="bold", size = 15),
      legend.text = element_text(face="bold", size = 10)
    ) +
    ggtitle(title)
  return(output.plot)
}