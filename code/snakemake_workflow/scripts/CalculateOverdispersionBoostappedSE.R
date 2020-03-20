setwd("/Users/benfair/Documents/GiladLiProjects/Repos/Comparative_eQTL/code/snakemake_workflow/")

library(tidyverse)
library(knitr)
library("edgeR")
source("../CustomFunctions.R")


GetLoessResidual <- function(x, y){
  loess.fit <- loess(y ~ x, degree=1)
  #for our purposes, the fit with the default degree=2 looks like unnatural fit, especially near the edges
  loess.residual <- y - predict(loess.fit, x)
  return(loess.residual)
}

### Pick samples to drop

CountTableChimpFile <- '../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz'
CountTableChimp <- read.table(gzfile(CountTableChimpFile), header=T, check.names=FALSE, skip=1)
colnames(CountTableChimp) <- paste0("C.", colnames(CountTableChimp))

CountTableHumanFile <- '../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz'
CountTableHuman <- read.table(gzfile(CountTableHumanFile), header=T, check.names=FALSE, skip=1)
colnames(CountTableHuman) <- paste0("H.", colnames(CountTableHuman))

CombinedTable <- inner_join(CountTableChimp[,c(1,7:length(CountTableChimp))], CountTableHuman[,c(1,7:length(CountTableHuman))], by=c("C.Geneid"="H.Geneid")) %>%
  column_to_rownames("C.Geneid") %>% as.matrix()

HumanSampleNames3 <- read.delim("../../output/PowerAnalysisCountTable.Human.SampleList.txt", col.names=c("SRA", "GTExName"), header=F, stringsAsFactors = F, sep=" ")

HumanSamplesToDropOutOfCluster <- c("SRR613186", "SRR598509", "SRR1478149", "SRR603918", "SRR1507229")
HumanSamplesToDropLowReads <- CombinedTable %>% colSums() %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(SampleNoSpecies = gsub("^[HC].", "", Sample)) %>%
  inner_join(HumanSampleNames3, by=c("SampleNoSpecies"="SRA")) %>%
  mutate(Species=substr(Sample, 1,1)) %>%
  filter(Species=="H") %>%
  filter(!SampleNoSpecies %in% HumanSamplesToDropOutOfCluster) %>%
  top_n(-5, wt=as.numeric(`.`)) %>% pull(SampleNoSpecies)
HumanSamplesToDrop <- c(HumanSamplesToDropOutOfCluster, HumanSamplesToDropLowReads)

ChimpSamplesToDrop <- c()
