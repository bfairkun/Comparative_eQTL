library(tidyverse)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
CountTableIn <- args[1]
SampleListIn <- args[2]
CountTableOut <- args[3]

# setwd("/project2/gilad/bjf79_project1/projects/Comparative_eQTL/code/snakemake_workflow/")
# CountTableIn <- "OverdispersionAnalysis/FullCountTable.txt.gz"
# SampleListIn <- "../../data/OverdispersionGTExAnalysisSampleList.txt"
# CountTableOut <- "scratch/CountTableSubset.txt"


SampleList <- read.table(SampleListIn, header=T, sep='\t', stringsAsFactors = F)
CountTable <-fread(paste('zcat', CountTableIn), sep='\t', header=T, skip=2)

CountTable %>% dplyr::select(Name, SampleList$SAMPID) %>%
  write.table(CountTableOut, row.names = F, sep='\t', quote=F)
