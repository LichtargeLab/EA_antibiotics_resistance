# EA_antibiotics_resistance

This repository contains all the code to analyze the driver genes in ciprofloxacin and colistin resistant E. coli using Evolutionary Action (EA). This method can be generalized to detect genes under other source of selections in evolved E. coli strains. To perform these analyses or query EA scores for mutations in MG1655, please visit our [R shiny app](http://bioheat.lichtargelab.org).

#### Websites

[Lichtarge Lab](http://lichtargelab.org)

[R shiny app to analyze driver genes in E. coli](http://bioheat.lichtargelab.org)

##### Related References

[Katsonis P, Lichtarge O., A formal perturbation equation between genotype and phenotype determines the evolutionary action of protein coding variations on fitness, Genome Research, 2014 Sep 12. pii: gr.176214.114.](https://pubmed.ncbi.nlm.nih.gov/25217195/)

### Dependencies

```
install.packages(c("tidyverse", "furrr", "RobustRankAggreg", "openxlsx", "scales", "cowplot", "ggh4x"))
devtools::install_github("rensa/stickylabeller")

```
