## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----load-libraries, message=F, warning=F--------------------------------
library(tidyverse)
library(knitr)


## ------------------------------------------------------------------------
# knitr::purl(input="../analysis/20190327_MakeFamPhenotypeFile.Rmd", output="../analysis/20190327_MakeFamPhenotypeFile.R")


## ----Set-filepaths-------------------------------------------------------
# Use command line input to specify input and output if this is the Rscript version of this file (as opposed to Rmarkdown).
if(commandArgs()[4] == "--file=../../analysis/20190327_MakeFamPhenotypeFile.R"){
  args <- commandArgs(trailingOnly = T)
  CountFilepath <- args[1]
  EmptyFamFilepath <- args[2]
  PhenotypeOutFilepath <- args[3]
  PhenotypeListOutFilepath <- args[4]
  GenesBedFile <- args[5]
} else {
  CountFilepath <- '../output/CountTable.tpm.txt.gz'
  EmptyFamFilepath <- '../output/ForAssociationTesting.temp.fam'
  GenesBedFile <- '../data/cDNA.all.chromosomal.bed'
}


## ----make-tidy-data, warning=F-------------------------------------------
CountTable <- read.table(gzfile(CountFilepath), header=T, check.names=FALSE, row.names = 1)

# dimensions of count table
CountTable %>% dim()

kable(head(CountTable))

EmptyFamFile <- read.table(EmptyFamFilepath, col.names=c("FID", "IID", "Father", "Mother", "SX", "Pheno"), stringsAsFactors = F) %>%
  select(-Pheno)

# Will use GeneRegions to filter out non-autosomal genes
GeneChromosomes <- read.table(GenesBedFile, col.names=c("chromosome", "start", "stop", "gene", "score", "strand"), stringsAsFactors = F) %>%
  select(gene, chromosome)
kable(head(GeneChromosomes))


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
  merge(GeneChromosomes, by="gene", all.x=T) %>%
  filter(!chromosome %in% c("X", "Y", "MT")) %>%
  select(-chromosome) %>%
  column_to_rownames('gene') %>%
  log() %>%
  t()

row.names(PhenotypesToOutput) <- colnames(CountTable)


Output.df <- EmptyFamFile %>%
  merge(PhenotypesToOutput, all.x=T, by.x="IID", by.y=0) %>% as.tibble()
Output.df

GeneList <- data.frame(GeneList=colnames(Output.df)[-1:-5])


## ----write-table-if-script-----------------------------------------------

if(commandArgs()[4] == "--file=../../analysis/20190327_MakeFamPhenotypeFile.R"){
  write.table(Output.df, col.names = F, sep='\t', file=PhenotypeOutFilepath, row.names=F, quote=F)
  write.table(GeneList, col.names = F, sep='\t', file=PhenotypeListOutFilepath, row.names=F, quote=F)
}


