library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
file_in <- args[1]
file_out <- args[2]

# setwd("/project/gilad/bjf79/projects/Comparative_eQTL/code/snakemake_workflow/")
# file_in <- "../../data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz"
# file_out <- "scratch/GtexTissueMatrixFilteredTest.tsv"

GTEx.Expression <- read.table(file_in, sep='\t', skip=2, header=T,check.names=FALSE)

#Get list of brain subtypes to exclude (let's just include brain cerebellum for this analysis)
#...too many brain subtypes is redundant
BrainsSubtypesToExclude<-colnames(GTEx.Expression) %>% grep("Brain", ., value=T) %>% grep("Cerebellum", ., value=T, invert=T)

GTEx.Expression %>%
  mutate(human_id=gsub("\\.\\d+$", "", gene_id, perl=T)) %>%
  select(-c(gene_id, BrainsSubtypesToExclude, Description)) %>%
  select(human_id, everything()) %>%
  write.table(file=file_out, quote=F, sep='\t', row.names = F)
