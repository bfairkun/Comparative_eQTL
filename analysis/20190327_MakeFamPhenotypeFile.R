## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----load-libraries, message=F, warning=F--------------------------------
library(tidyverse)
library(knitr)
library(matrixStats)



## ------------------------------------------------------------------------
# knitr::purl(input="../analysis/20190327_MakeFamPhenotypeFile.Rmd", output="../analysis/20190327_MakeFamPhenotypeFile.R")


## ----Set-filepaths-------------------------------------------------------
# Use command line input to specify input and output if this is the Rscript version of this file (as opposed to Rmarkdown).
if(commandArgs()[4] == "--file=../../analysis/20190327_MakeFamPhenotypeFile.R"){
  args <- commandArgs(trailingOnly = T)
  CountFilepath <- args[1]
  EmptyFamFilepath <- args[2]
  PhenotypeOutFilepath <- args[3]
  GenesBedFile <- args[4]
  TranscriptToGenesFile <- args[5] #Two column mapping of gene_id(col1) and transcript_id(col2)
} else {
  CountFilepath <- '../output/CountTable.tpm.txt.gz'
  EmptyFamFilepath <- '../output/ForAssociationTesting.temp.fam'
  GenesBedFile <- '../output/Genes.bed'
  TranscriptToGenesFile <- '../data/Biomart_export.Pan_Tro_3.geneids.txt'
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

# ensembl transcript and gene-ids
TranscriptsToGenes <- read.table(TranscriptToGenesFile, sep='\t', header=T, stringsAsFactors = F) %>% select(Gene.stable.ID, Transcript.stable.ID)
kable(head(TranscriptsToGenes))


## ------------------------------------------------------------------------
#How many genes are there
length(unique(TranscriptsToGenes$Gene.stable.ID))

#filter count table for autosomal genes, and sum transcript TPM to gene level TPM.
MyDf <- CountTable %>%
  rownames_to_column('gene') %>%
  mutate(Transcript.stable.ID=sub(".\\d+$", "", gene)) %>%
  select(-gene) %>%
  left_join(TranscriptsToGenes, by="Transcript.stable.ID") %>%
  select(-Transcript.stable.ID) %>%
  aggregate(.~Gene.stable.ID, ., sum) %>%
  left_join(GeneChromosomes, by=c("Gene.stable.ID"="gene")) %>%
  filter(!chromosome %in% c("X", "Y", "MT")) %>%
  select(-chromosome)

  
Median <- MyDf %>%
  column_to_rownames('Gene.stable.ID') %>%
  mutate(Median = rowMedians(as.matrix(.))) %>% pull(Median)

# histogram of median expression [ log10(TPM) ] of genes 
hist(log10(Median))

GeneSet1 <- MyDf %>%
  mutate(Median=Median) %>%
  filter(Median>0.1) %>%
  pull(Gene.stable.ID)

GeneSet2 <- MyDf %>%
  filter_if(is.numeric, all_vars(.>0)) %>%
  pull(Gene.stable.ID)

#Number genes left after filter method1
length(GeneSet1)

#Number genes left after filter method2
length(GeneSet2)

#Number genes left after intersection of both methods
length(intersect(GeneSet1, GeneSet2))



## ------------------------------------------------------------------------

# Output is log10(TPM)
PhenotypesToOutput <- MyDf %>%
  filter_if(is.numeric, all_vars(.>0)) %>%
  column_to_rownames('Gene.stable.ID') %>%
  log10()

# row.names(PhenotypesToOutput) <- colnames(CountTable)
# 
# 
# Output.df <- EmptyFamFile %>%
#   merge(PhenotypesToOutput, all.x=T, by.x="IID", by.y=0) %>% as.tibble()
# Output.df
# 
# GeneList <- data.frame(GeneList=colnames(Output.df)[-1:-5])


## ----write-table-if-script-----------------------------------------------

if(commandArgs()[4] == "--file=../../analysis/20190327_MakeFamPhenotypeFile.R"){
  write.table(PhenotypesToOutput, col.names = NA, sep='\t', file=PhenotypeOutFilepath, row.names=T, quote=F)
}


