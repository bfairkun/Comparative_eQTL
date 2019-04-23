## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----load-libraries, message=F, warning=F--------------------------------
library(plyr)
library(tidyverse)
library(knitr)
library(readxl)


## ------------------------------------------------------------------------
# knitr::purl(input="../analysis/20190327_MakeCovariateFiles.Rmd", output="../analysis/20190327_MakeCovariateFiles.R")


## ----Set-filepaths-------------------------------------------------------
# Use command line input to specify input and output if this is the Rscript version of this file (as opposed to Rmarkdown).
if(commandArgs()[4] == "--file=../../analysis/20190327_MakeCovariateFiles.R"){
  args <- commandArgs(trailingOnly = T)
  EmptyFamFilepath <- args[1]
  MetadataExcel <- args[2]
  CountFilepath <- args[3]
  NumRS.PCs <- args[4]
  GenotypePCs <- args[5]
  NumGenotype.PCs <- args[6]
  OutputFile <- args[7]
  print (args)

} else {
  CountFilepath <- '../output/CountTable.tpm.txt.gz'
  EmptyFamFilepath <- '../output/ForAssociationTesting.temp.fam'
  MetadataExcel <- '../data/Metadata.xlsx'
  GenotypePCs <- '../output/PopulationStructure/pca.eigenvec'
  NumGenotype.PCs <- 3
  NumRS.PCs <- 3
}


## ----make-tidy-data, warning=F-------------------------------------------
CountTable <- read.table(gzfile(CountFilepath), header=T, check.names=FALSE, row.names = 1)

# dimensions of count table
CountTable %>% dim()

kable(head(CountTable))

EmptyFamFile <- read.table(EmptyFamFilepath, col.names=c("FID", "IID", "Father", "Mother", "SX", "Pheno"), stringsAsFactors = F) %>%
  select(-Pheno)

GenotypePCs.df <- read.table(GenotypePCs, header=T) %>% 
  select(-c("FID")) %>%
  rename_all(paste0, ".GT")
kable(head(GenotypePCs.df))

OtherMetadata <- as.data.frame(read_excel(MetadataExcel))
kable(head(OtherMetadata))


## ----write-covariates----------------------------------------------------
# RNA-seq PCA
pca_results <- CountTable %>%
  +0.1 %>%
  mutate(sumVar = rowSums(.)) %>%
  arrange(desc(sumVar)) %>%
  head(2500) %>%
  select(-sumVar) %>%
  log() %>%
  t() %>%
  prcomp(center=T, scale. = T)

RNASeqPCs.df <- pca_results$x %>% 
  as.data.frame %>% 
  rownames_to_column() %>%
  rename_all(paste0, ".RS")

Covariates <- EmptyFamFile %>%
  select(IID) %>%
  left_join(RNASeqPCs.df[0:NumRS.PCs+1], by=c("IID" = "rowname.RS")) %>%
  left_join(GenotypePCs.df[0:NumGenotype.PCs+1], by=c("IID" = "IID.GT")) %>%
  left_join(
    OtherMetadata %>% select(Individual.ID, SX),
    by = c("IID" = "Individual.ID")
  ) %>%
  mutate(SX = plyr::mapvalues(SX, from=c("M", "F"), to=c(0,1))) %>%
  t()
  
kable(head(Covariates))


## ----write-table-if-script-----------------------------------------------
if(commandArgs()[4] == "--file=../../analysis/20190327_MakeCovariateFiles.R"){
  write.table(Covariates, col.names = F, sep='\t', file=OutputFile, row.names=T, quote=F)
}

