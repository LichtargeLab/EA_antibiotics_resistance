# Method adapted by Chen Wang from Carr, I. M., J. I. Robinson, R. Dimitriou, 
# A. F. Markham, A. W. Morgan, and D. T. Bonthron. 2009. “Inferring Relative 
# Proportions of DNA Variants from Sequencing Electropherograms.” Bioinformatics 
# 25 (24): 3244–50. https://doi.org/10.1093/bioinformatics/btp583.
# This script computes the allele frequency at a give position using Sanger
# sequencing data. Please refer to the original publication for detail.

# Check and install packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if ("sangerseqR" %in% installed.packages() == FALSE) {
  BiocManager::install("sangerseqR")
}
if ("dplyr" %in% installed.packages() == FALSE) {
  install.packages("dplyr")
}

# Load packages
library(sangerseqR)
library(dplyr)

# Used to locate the mutation site
dna_up <- DNAString("ATTCTCGGTC") # Xnt upstream of the target mutation site. 
# Make sure to use the complement stand when using reverse primer in sequencing.

wt <- "C" # wt version
mut <- "T" # mut version
mypath <- "AF_from_SangerSeq_example/" # The folder where the electropherogram are saved.
# Two reference electropherograms that have pure wt and pure mut should be included in this folder.
# A csv file will be generated in the same folder.

# Uses -10 ~ -6 peaks as references
test.files <- list.files(path = mypath, pattern = "*.ab1", full.names = FALSE)
wt.peak.adj <- rep(0,length(test.files))
mut.peak.adj <- rep(0,length(test.files))
posi <- rep(0,length(test.files))
note <- rep("",length(test.files))
for (i in 1:length(test.files)) {
  test <- readsangerseq(paste0(mypath, test.files[i]))
  raw.data <- traceMatrix(test)[peakPosMatrix(test)[,1],]
  colnames(raw.data) <- c("A", "C", "G", "T")
  firstpeaks <- rep(0, nrow(raw.data))
  for (j in 1:nrow(raw.data)) {
    firstpeaks [j] <- max(raw.data[j,])
  }
  
  dna.start <- matchPattern(dna_up, primarySeq(test))@ranges@start
  if (length(dna.start)==0) {
    wt.peak.adj[i] <- NA
    mut.peak.adj[i] <- NA
    posi[i] <- NA
    next
  }
  mut.pos <- dna.start + matchPattern(dna_up, primarySeq(test))@ranges@width
  ref.start <- mut.pos - 10
  ref.end <- ref.start + 4
  
  # Within each sample, adjust the peak height of position of 
  # interest to the peak for positions -10 to -6
  ref.peak <- mean(firstpeaks[ref.start:ref.end])
  wt.peak.adj[i] <- as.numeric(raw.data[mut.pos,wt]/ref.peak)
  mut.peak.adj[i] <- as.numeric(raw.data[mut.pos,mut]/ref.peak)
  posi[i] <- mut.pos
  
  # Check the main call for target position in the electropherogram.
  # Throws a warning in the note column if neither wt nor mut is the main call.
  # This might indicate an mislabel of the wt/ref base or poor sequencing quality.
  main_call <- names(which.max(raw.data[mut.pos,]))
  if (! main_call %in% c(wt, mut)) {
    note[i] <- paste0("The main call is ", main_call, ", not wt or mut.")
  }
}
work.sheet <- tibble(test.files, wt.peak.adj, mut.peak.adj, posi, note)
output <- work.sheet %>% 
  filter(is.na(posi)==FALSE) %>%
  mutate(wt.ref = case_when(wt.peak.adj == min(wt.peak.adj) ~ "min",
                            wt.peak.adj == max(wt.peak.adj) ~ "max",
                            .default = "")) %>%
  mutate(mut.ref = case_when(mut.peak.adj == min(mut.peak.adj) ~ "min",
                            mut.peak.adj == max(mut.peak.adj) ~ "max",
                            .default = "")) %>%
  mutate(wt.peak.adj=wt.peak.adj-min(wt.peak.adj)) %>%
  mutate(mut.peak.adj=mut.peak.adj-min(mut.peak.adj)) %>%
  mutate(wt.peak.adj=wt.peak.adj/max(wt.peak.adj)) %>%
  mutate(mut.peak.adj=mut.peak.adj/max(mut.peak.adj)) %>%
  mutate(AF=mut.peak.adj/(mut.peak.adj+wt.peak.adj)) %>%
  relocate(test.files, wt.peak.adj, mut.peak.adj, AF, posi, wt.ref, mut.ref, note) %>%
  as.data.frame() 

# Usually the pure wt sample is picked as wt max reference and mut min reference
# Similarly the pure mut sample is picked as mut max reference and wt min reference

output.NA <- work.sheet %>% 
  filter(is.na(posi)==TRUE)
output <- bind_rows(output, output.NA)
write.csv(output, file = paste0(mypath, "AF_analysis.csv"), row.names = FALSE)
