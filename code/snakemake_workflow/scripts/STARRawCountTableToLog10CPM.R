library(tidyverse)
library(matrixStats)


args <- commandArgs(trailingOnly = T)
CountFilepath <- args[1]
EmptyFamFilepath <- args[2]
PhenotypeOutFilepath <- args[3]
GenesBedFile <- args[4]

EmptyFamFile <- read.table(EmptyFamFilepath, col.names=c("FID", "IID", "Father", "Mother", "SX", "Pheno"), stringsAsFactors = F) %>%
  select(-Pheno)

## ----make-tidy-data, warning=F-------------------------------------------
CountTable <- read.table(gzfile(CountFilepath), header=T, check.names=FALSE, row.names = 1)

# dimensions of count table
CountTable %>% dim()

head(CountTable)

EmptyFamFile <- read.table(EmptyFamFilepath, col.names=c("FID", "IID", "Father", "Mother", "SX", "Pheno"), stringsAsFactors = F) %>%
  select(-Pheno)

# Will use GeneRegions to filter out non-autosomal genes
GeneChromosomes <- read.table(GenesBedFile, col.names=c("chromosome", "start", "stop", "gene", "score", "strand"), stringsAsFactors = F) %>%
  select(gene, chromosome)
head(GeneChromosomes)


## ------------------------------------------------------------------------

dim(CountTable)
#filter count table for autosomal genes, and sum transcript TPM to gene level TPM.
MyDf <- CountTable %>%
  rownames_to_column('gene') %>%
  left_join(GeneChromosomes, by=c("gene")) %>%
  filter(!chromosome %in% c("X", "Y", "MT")) %>%
  select(-chromosome) %>%
  column_to_rownames('gene')

head(MyDf)

Q <- rowQuantiles(as.matrix(MyDf), probs=0.2)
# Now some filtering:
GeneSetFilter1 <- names(which(Q>10))

CountTableFiltered <- MyDf %>% as.data.frame() %>%
  rownames_to_column('gene') %>%
  filter(gene %in% GeneSetFilter1) %>%
  column_to_rownames('gene')

# table of log10CPM
log10CPM_table <- log10((CountTableFiltered + 1) / colSums(CountTableFiltered))
dim(log10CPM_table)

write.table(log10CPM_table, col.names = NA, sep='\t', file=PhenotypeOutFilepath, row.names=T, quote=F)
