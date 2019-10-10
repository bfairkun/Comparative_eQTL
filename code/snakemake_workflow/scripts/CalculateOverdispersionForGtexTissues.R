library(tidyverse)
library(data.table)
source("../CustomFunctions.R")

args = commandArgs(trailingOnly=TRUE)
CountTableIn <- args[1]
SampleListIn <- args[2]
GeneLengthsIn <- args[3]
MuTableOut <- args[4]
OverdispersionTableOut <- args[5]


# setwd("/project2/gilad/bjf79_project1/projects/Comparative_eQTL/code/snakemake_workflow/")
# CountTableIn <- "OverdispersionAnalysis/TableSubset.txt"
# SampleListIn <- "../../data/OverdispersionGTExAnalysisSampleList.txt"
# GeneLengthsIn <- "OverdispersionAnalysis/GTEx.genelengths.txt"
# MuTableOut <- args[3]
# OverdispersionTableOut <- args[4]


SampleList <- read.table(SampleListIn, header=T, sep='\t', stringsAsFactors = F)
CountTable <-fread(CountTableIn, sep='\t', header=T)
GeneLengths <- read.table(GeneLengthsIn, header=T, sep='\t', stringsAsFactors = F)


GeneLengthVector <- CountTable %>% dplyr::select("Name") %>%
  left_join(GeneLengths, by=c("Name"="gene")) %>%
  pull(mean)

TissueList <- SampleList %>% distinct(SMTSD) %>% pull(SMTSD)


MuMatrix <- matrix(nrow=length(GeneLengthVector), ncol=length(TissueList))
colnames(MuMatrix) <- TissueList
rownames(MuMatrix) <- CountTable$Name

OverdispersionMatrix <- MuMatrix

for (i in seq_along(TissueList)){
  print(paste("Analyzing", TissueList[i]))
  SamplesToSubset <- SampleList %>% filter(SMTSD == TissueList[i]) %>% pull(SAMPID)
  CountTableSubset <- CountTable %>%
    dplyr::select(Name, SamplesToSubset) %>%
    column_to_rownames("Name")
  ParameterEstimates <- GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTableSubset, GeneLengthVector)
  MuMatrix[,i] <- ParameterEstimates$mu
  OverdispersionMatrix[,i] <- ParameterEstimates$overdispersion
}

MuMatrix %>% as.data.frame() %>%
  write.table(MuTableOut, sep='\t', quote = F, col.names = NA)

# MuMatrix %>% as.data.frame() %>%
#   write.table("../../output/GTEx.Tissues.Mu.Matrix.txt", sep='\t', quote = F, col.names = NA)


OverdispersionMatrix %>% as.data.frame() %>%
  write.table(OverdispersionTableOut, sep='\t', quote = F, col.names = NA)

# OverdispersionMatrix %>% as.data.frame() %>%
#   write.table("../../output/GTEx.Tissues.Overdispersion.Matrix.txt", sep='\t', quote = F, col.names = NA)
# heatmap.2(B, trace="none", dendrogram=c("col"), labCol=F, cexRow = 0.7)
# heatmap.2(A, trace="none", dendrogram=c("col"), labCol=F, cexRow = 0.7)

