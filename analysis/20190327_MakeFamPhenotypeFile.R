## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load-libraries, message=F, warning=F--------------------------------
library(tidyverse)
library(knitr)

## ----Set-filepaths-------------------------------------------------------
# Use command line input to specify input and output if this is the Rscript version of this file (as opposed to Rmarkdown).
if(commandArgs()[4] == "--file=../../analysis/20190327_MakeFamPhenotypeFile.R"){
  args <- commandArgs(trailingOnly = T)
  CountFilepath <- args[1]
  EmptyFamFilepath <- args[2]
  PhenotypeOutFilepath <- args[3]
} else {
  CountFilepath <- '../output/CountTable.tpm.txt.gz'
  EmptyFamFilepath <- '../output/ForAssociationTesting.temp.fam'
}

## ----make-tidy-data, warning=F-------------------------------------------
CountTable <- read.table(gzfile(CountFilepath), header=T, check.names=FALSE, row.names = 1)


# dimensions of count table
CountTable %>% dim()

kable(head(CountTable))

EmptyFamFile <- read.table(EmptyFamFilepath, col.names=c("FID", "IID", "Father", "Mother", "SX", "Pheno"), stringsAsFactors = F) %>%
  select(-Pheno)

## ------------------------------------------------------------------------
GeneSet1 <- CountTable %>%
  rownames_to_column('gene') %>%
  mutate(Mean = rowMeans(select(., -gene))) %>%
  filter(Mean>1) %>%
  pull(gene)

GeneSet2 <- CountTable %>%
  rownames_to_column('gene') %>%
  filter_if(is.numeric, all_vars(.>0)) %>%
  pull(gene)

#Number genes left after filter method1
length(GeneSet1)

#Number genes left after filter method2
length(GeneSet2)

#Number genes left after intersection of both methods
length(intersect(GeneSet1, GeneSet2))

## ------------------------------------------------------------------------
PhenotypesToOutput <- CountTable %>%
  rownames_to_column('gene') %>%
  filter_if(is.numeric, all_vars(.>0)) %>%
  column_to_rownames('gene') %>%
  log() %>%
  t()

row.names(PhenotypesToOutput) <- colnames(CountTable)

Output.df <- EmptyFamFile %>%
  merge(PhenotypesToOutput, all.x=T, by.x="IID", by.y=0) %>% as.tibble()
Output.df

## ----write-table-if-script-----------------------------------------------

if(commandArgs()[4] == "--file=../../analysis/20190327_MakeFamPhenotypeFile.R"){
  write.table(Output.df, col.names = F, sep='\t', output=PhenotypeOutFilepath)
}


