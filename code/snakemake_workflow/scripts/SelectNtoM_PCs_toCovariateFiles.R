library(tidyverse)

args <- commandArgs(trailingOnly = T)
PCsFromLeafcutter <- args[1]
EmptyFamFilepath <- args[2]
PhenotypeOutBaseFilepath <- args[3]
N <- args[4]
M <- args[5]

EmptyFamFile <- read.table(EmptyFamFilepath, col.names=c("FID", "IID", "Father", "Mother", "SX", "Pheno"), stringsAsFactors = F, sep=" ") %>%
  select(-Pheno)

Covariates <- read.table(PCsFromLeafcutter, header=T, check.names = F)

for (i in N:M){
  Covariates %>%
    select(id, EmptyFamFile$IID) %>%
    filter(id %in% N:i) %>%
    write.table(file=paste0(PhenotypeOutBaseFilepath, i, ".covariates.txt"), quote=F, row.names=F, sep='\t')
}
