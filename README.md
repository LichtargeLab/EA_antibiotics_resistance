# EA_antibiotics_resistance

This repository contains all the code to analyze the driver genes in ciprofloxacin and colistin resistant E. coli using Evolutionary Action (EA). This method can be generalized to detect genes under other source of selections in evolved E. coli strains. To perform these analyses or query EA scores for mutations in MG1655, please visit our [R shiny app](http://bioheat.lichtargelab.org).

Raw sequencing files for the adaptive evolution experiments can be accessed [here](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA543834).

Use [src/AF_from_SangerSeq.R](https://github.com/LichtargeLab/EA_antibiotics_resistance/blob/master/src/AF_from_SangerSeq.R) to compute allele frequencies from electropherograms.

### Websites

[Lichtarge Lab](http://lichtargelab.org)

[R shiny app to analyze driver genes in E. coli](http://bioheat.lichtargelab.org)

### Related References

[Marciano DC, Wang C, Hsu TK, Bourquard T, Atri B, Nehring RB, Abel NS, Bowling EA, Chen TJ, Lurie PD, Katsonis P, Rosenberg SM, Herman C, Lichtarge O. Evolutionary action of mutations reveals antimicrobial resistance genes in Escherichia coli. Nat Commun. 2022 Jun 9;13(1):3189.](https://www.nature.com/articles/s41467-022-30889-1)

[Katsonis P, Lichtarge O. A formal perturbation equation between genotype and phenotype determines the Evolutionary Action of protein-coding variations on fitness. Genome Res. 2014 Dec;24(12):2050-8.](https://genome.cshlp.org/content/24/12/2050.long)

[Carr IM, Robinson JI, Dimitriou R, Markham AF, Morgan AW, Bonthron DT. Inferring relative proportions of DNA variants from sequencing electropherograms. Bioinformatics. 2009 Dec 15;25(24):3244-50.](https://academic.oup.com/bioinformatics/article/25/24/3244/235657?login=false)

### Dependencies

```
install.packages(c("tidyverse", "furrr", "RobustRankAggreg", "openxlsx", "scales", "cowplot", "ggh4x"))
devtools::install_github("rensa/stickylabeller")

```

### Notes

The code is tested on macOS Catalina. The EA analyses took several minutes to run on a 2017 MacBook Air for both ALE and natural isolate datasets. The downsampling analyses took days to run on a computing cluster using 20 computing threads.