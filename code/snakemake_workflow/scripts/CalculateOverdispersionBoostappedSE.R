# setwd("/Users/benfair/Documents/GiladLiProjects/Repos/Comparative_eQTL/code/snakemake_workflow/")

library(tidyverse)
library("edgeR")
source("../CustomFunctions.R")


GetLoessResidual <- function(x, y){
  loess.fit <- loess(y ~ x, degree=1)
  #for our purposes, the fit with the default degree=2 looks like unnatural fit, especially near the edges
  loess.residual <- y - predict(loess.fit, x)
  return(loess.residual)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

Resampled.GetDispersion <- function(CountMatrix, GeneLengthsMatrix, n=NULL, N_genes=NULL, seed=0){
  set.seed(seed)
  if(is.null(N_genes)){
    NumGenes <- nrow(CountMatrix)
  }
  else{
    NumGenes<-N_genes
  }
  if(is.null(n)){
    NumInds <- ncol(CountMatrix)
  }
  else{
    NumInds<-n
  }
  ColsSample <- sample(x=ncol(CountMatrix), size=NumInds, replace=T)
  Estimates <- GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountMatrix[1:NumGenes,ColsSample], GeneLengthsMatrix[1:NumGenes,ColsSample])
  Dispersion <- GetLoessResidual(Estimates$mu, log(Estimates$overdispersion))

  return(setNames(Dispersion, rownames(CountMatrix[1:NumGenes,])))
}

### Global vars
CountTableChimpFile <- '../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz'
CountTableHumanFile <- '../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz'
OutputDE <- '../../output/Final/TableS1.tab'
OutputDispersionEstimatesChimp <- 'OverdispersionBootsrappedIterations/Chimp.1.2.tab'
OutputDispersionEstimatesHuman <- 'OverdispersionBootsrappedIterations/Human.1.2.tab'
Npermutations <- 5
InitialSeed <- 1

args = commandArgs(trailingOnly=TRUE)
CountTableChimpFile <- args[1]
CountTableHumanFile <- args[2]
OutputDE <- args[3]
OutputDispersionEstimatesChimp <- args[4]
OutputDispersionEstimatesHuman <- args[5]
Npermutations <- args[6]
InitialSeed <- args[7]

### Pick samples to drop

HumanSamplesToDrop <- c("SRR613186",  "SRR598509",  "SRR1478149", "SRR603918",  "SRR1507229", "SRR1478900", "SRR1477015", "SRR601986",  "SRR614996",  "SRR1474730")
ChimpSamplesToDrop <- c()

### Pick genes to analyze

DE.results <- read.delim(OutputDE, sep='\t', stringsAsFactors = F)
GeneListForOverdispersionCalculation <- DE.results$Ensembl_geneID

NumRowsToAnalyze=length(GeneListForOverdispersionCalculation)
NumRowsToAnalyze=100

### Estimate overdispersion after random sampling with replacement of columns
CountTables <- GetCountTables(CountTableChimpFile,
                              CountTableHumanFile,
                              0, GeneListForOverdispersionCalculation, ChimpSampleDrop=ChimpSamplesToDrop, HumanSampleDrop = HumanSamplesToDrop)


HumanResultsMatrix <-matrix(data=NA, nrow=NumRowsToAnalyze, ncol=Npermutations)
ChimpResultsMatrix <-matrix(data=NA, nrow=NumRowsToAnalyze, ncol=Npermutations)

for (i in 1:Npermutations){
  HumanResultsMatrix[,i] <- Resampled.GetDispersion(CountTables$Human$Counts,
                                                    rep.col(CountTables$Human$GeneLengths, ncol(CountTables$Human$Counts)),
                                                    seed=i,
                                                    n=39,
                                                    N_genes=NumRowsToAnalyze)
  ChimpResultsMatrix[,i] <- Resampled.GetDispersion(CountTables$Chimp$Counts,
                                                    rep.col(CountTables$Chimp$GeneLengths, ncol(CountTables$Chimp$Counts)),
                                                    seed=i,
                                                    n=39,
                                                    N_genes=NumRowsToAnalyze)
}

row.names(ChimpResultsMatrix) <- rownames(CountTables$Chimp$Counts[1:NumRowsToAnalyze,])
row.names(HumanResultsMatrix) <- rownames(CountTables$Human$Counts[1:NumRowsToAnalyze,])

write.table(t(ChimpResultsMatrix), OutputDispersionEstimatesChimp, quote=F, sep='\t', col.names = T, row.names=F)
write.table(t(HumanResultsMatrix), OutputDispersionEstimatesHuman, quote=F, sep='\t', col.names = T, row.names=F)


