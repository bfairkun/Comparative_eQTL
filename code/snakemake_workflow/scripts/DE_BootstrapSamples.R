library(tidyverse)
library("edgeR")
source("../CustomFunctions.R")

# InitialSeed <- 0
# CountTableChimpSubreadFile <- '../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz'
# CountTableHumanSubreadFile <- '../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz'
# OutputDepthName <- "test"
# OutputFn <- "scratch/Test.DE.txt.gz"


args = commandArgs(trailingOnly=TRUE)
InitialSeed <- as.integer(args[1])
CountTableChimpSubreadFile <- args[2]
CountTableHumanSubreadFile <- args[3]
OutputDepthName <- args[4]
OutputFn <-args[5]


CountTableChimpFile <- '../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz'
CountTableHumanFile <- '../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz'

SampleSizesToTest <- c(2,4,8,12,16,24,32,39)

DropFile <- read.delim("../../data/DE_SamplesToDrop.txt", sep='\t', col.names = c("Sample", "Species"), stringsAsFactors = F)

HumanSamplesToDrop <- DropFile %>% filter(Species=="Human") %>% pull(Sample)
ChimpSamplesToDrop <- DropFile %>% filter(Species=="Chimp") %>% pull(Sample)

CountTableChimp <- read.table(gzfile(CountTableChimpFile), header=T, check.names=FALSE, skip=1)
colnames(CountTableChimp) <- paste0("C.", colnames(CountTableChimp))

CountTableHuman <- read.table(gzfile(CountTableHumanFile), header=T, check.names=FALSE, skip=1)
colnames(CountTableHuman) <- paste0("H.", colnames(CountTableHuman))

CombinedTable <- inner_join(CountTableChimp[,c(1,7:length(CountTableChimp))], CountTableHuman[,c(1,7:length(CountTableHuman))], by=c("C.Geneid"="H.Geneid")) %>%
  column_to_rownames("C.Geneid") %>%
  dplyr::select(-c(paste0("H.", HumanSamplesToDrop))) %>%
  as.matrix()

cpm <- cpm(CombinedTable, log=TRUE, prior.count=0.5)
d0 <- DGEList(CombinedTable)

#Calculate normalization factors
d0 <- calcNormFactors(d0)
cutoff <- 6
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d)

FC.NullInterval <- log2(1.0)

set.seed(InitialSeed)
TableOut <- data.frame()
for (i in SampleSizesToTest){
  print(i)
  True.efit <- DE.Subsampled(CountTableChimpSubreadFile, CountTableHumanSubreadFile, SubsampleSize = i, FC.NullInterval = FC.NullInterval, drop=drop, ChimpSampleDrop = NULL, HumanSampleDrop = HumanSamplesToDrop, Replacement = T)
  TableOut <- toptable(True.efit, n=Inf, sort.by="none") %>%
    rownames_to_column(var="Ensembl.ID") %>%
    dplyr::select(Ensembl.ID, logFC, P.Value, adj.P.Val) %>%
    mutate(SampleSize=i,
           Seed=InitialSeed,
           DepthFactor=OutputDepthName,
           logFC=signif(logFC,4),
           P.Value=signif(P.Value,2),
           adj.P.Val=signif(adj.P.Val,2)) %>%
    bind_rows(TableOut)
}

write_delim(TableOut, path=OutputFn, col_names=F, delim='\t')

