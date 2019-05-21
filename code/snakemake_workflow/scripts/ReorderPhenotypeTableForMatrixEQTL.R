library(tidyverse)

args <- commandArgs(trailingOnly = T)
PhenotypeTableFilepath <- args[1]
EmptyFamFilepath <- args[2]
PhenotypeOutFilepath <- args[3]
EmptyFamFile <- read.table(EmptyFamFilepath, col.names=c("FID", "IID", "Father", "Mother", "SX", "Pheno"), stringsAsFactors = F) %>%
  select(-Pheno)

PhenotypeFile <- read.table(PhenotypeTableFilepath, header=T, check.names = F)
PhenotypeFile %>%
  select(ID, EmptyFamFile$IID) %>%
  write.table(file=PhenotypeOutFilepath, col.names = T, row.names = F, quote=F, sep='\t')

# write.table(Output.df, col.names = F, sep='\t', file=PhenotypeOutFilepath, row.names=F, quote=F)
# write.table(GeneList, col.names = F, sep='\t', file=GeneListOutFilepath, row.names=F, quote=F)
