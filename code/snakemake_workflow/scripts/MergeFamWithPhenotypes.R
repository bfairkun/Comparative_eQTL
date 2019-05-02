library(tidyverse)

args <- commandArgs(trailingOnly = T)
PhenotypeTableFilepath <- args[1]
EmptyFamFilepath <- args[2]
PhenotypeOutFilepath <- args[3]
GeneListOutFilepath <- args[4]

EmptyFamFile <- read.table(EmptyFamFilepath, col.names=c("FID", "IID", "Father", "Mother", "SX", "Pheno"), stringsAsFactors = F) %>%
  select(-Pheno)
EmptyFamFile$IID

PhenotypeFile <- read.table(PhenotypeTableFilepath, header=T, check.names = F, row.names = 1)
Transposed <- t(PhenotypeFile) %>% as.data.frame() %>% rownames_to_column("IID")

Output.df <- EmptyFamFile %>%
  left_join(Transposed, by="IID") %>% as_tibble()
Output.df

GeneList <- data.frame(GeneList=colnames(Output.df)[-1:-5])

write.table(Output.df, col.names = F, sep='\t', file=PhenotypeOutFilepath, row.names=F, quote=F)
write.table(GeneList, col.names = F, sep='\t', file=GeneListOutFilepath, row.names=F, quote=F)
