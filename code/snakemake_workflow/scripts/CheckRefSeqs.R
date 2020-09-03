library(tidyverse)

args =  commandArgs(trailingOnly=TRUE)
EnsemblFai_f = args[1]

EnsemblFai_f <- "/project2/gilad/bjf79/genomes/Pan_tro_3.0_Ensembl/Sequence/DNA/Pan_troglodytes.Pan_tro_3.0.dna_sm.chromosome.fa.fai"

EnsemblFai <- read.delim(EnsemblFai_f, header = F, skip=0, col.names = c("name", "length", "offset", "linebases", "linewidth"))
GCAReport <- read.delim("Misc/EnsemblRef/Ref.report.tab", skip=32, header=T,  stringsAsFactors = F) %>%
  filter(Sequence.Role=="assembled-molecule") %>%
  select(GenBank.Accn, Sequence.Length, name=X..Sequence.Name, Assigned.Molecule) %>%
  filter(Assigned.Molecule %in% EnsemblFai$name)

GCAReport %>% filter(Assigned.Molecule=="MT")

GCAReport %>% select(Assigned.Molecule, GenBank.Accn) %>%
  write_delim("MiscOutput/ChrConversionKey.tab", delim ='\t', col_names=F)


