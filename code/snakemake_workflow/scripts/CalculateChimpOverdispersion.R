setwd("/project2/gilad/bjf79_project1/projects/Comparative_eQTL/code/snakemake_workflow/")
# setwd("/Users/benfair/Documents/GiladLiProjects/Repos/Comparative_eQTL/code/snakemake_workflow/")

source("../CustomFunctions.R")
library(readxl)

GetLoessResidual <- function(x, y){
  loess.fit <- loess(y ~ x, degree=1)
  #for our purposes, the fit with the default degree=2 looks like unnatural fit, especially near the edges
  loess.residual <- y - predict(loess.fit, x)
  return(loess.residual)
}

HumanSamplesToDrop <- c(c("SRR1507229","SRR603918", "SRR1478149", "SRR598509", "SRR613186"), c("SRR1489693", "SRR598148", "59167", "SRR1478900", "SRR1474730", "61317"))
ChimpSamplesToDrop <- c("Little_R")
OtherMetadata <- as.data.frame(read_excel("../../data/Metadata.xlsx"))
VirusChallengedChimps <- OtherMetadata %>% filter(grepl("V+",Viral.status)) %>% pull(Individual.ID)

CountTableChimpFile <- '../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz'
CountTableHumanFile <- '../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz'

DE.results <- read.table("../../data/DE_genes.NoVirusChallangedInds.txt", sep='\t', header=T, stringsAsFactors = F)
GeneListForOverdispersionCalculation <- DE.results$gene

### estimate overdispersion
CountTables <- GetCountTables(CountTableChimpFile,
                              CountTableHumanFile,
                              0, GeneListForOverdispersionCalculation, ChimpSampleDrop=ChimpSamplesToDrop, HumanSampleDrop = HumanSamplesToDrop)

NumRowsToAnalyze=length(GeneListForOverdispersionCalculation)
# NumRowsToAnalyze=10


Chimp.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Chimp$Counts[1:NumRowsToAnalyze,], CountTables$Chimp$GeneLengths[1:NumRowsToAnalyze])
Human.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Human$Counts[1:NumRowsToAnalyze,], CountTables$Human$GeneLengths[1:NumRowsToAnalyze])

cbind(Chimp.NB.fit.parameters,
                apply(CountTables$Chimp$log2RPKM, 1, mean)[1:NumRowsToAnalyze],
                Human.NB.fit.parameters,
                apply(CountTables$Human$log2RPKM, 1, mean)[1:NumRowsToAnalyze]) %>%
  setNames(
    c("Chimp.Mean.Expression", "Chimp.Overdispersion", "Chimp.theta.se", "Chimp.Mean.Log2RPKM", "Human.Mean.Expression", "Human.Overdispersion", "Human.theta.se","Human.Mean.Log2RPKM")
  ) %>%
  rownames_to_column(var="gene") %>%
  mutate(
    Chimp.Residual = GetLoessResidual(Chimp.Mean.Expression, log(Chimp.Overdispersion)),
    Human.Residual = GetLoessResidual(Human.Mean.Expression, log(Human.Overdispersion))

  ) %>%
  mutate(
    Chimp.Dispersion = Chimp.Residual + median(log(Chimp.Overdispersion), na.rm = T),
    Human.Dispersion = Human.Residual + median(log(Human.Overdispersion), na.rm = T)
  ) %>%
  dplyr::select(gene, everything()) %>%
write.table("../../output/OverdispersionEstimatesFromChimp.txt", sep='\t', quote=F, row.names = F)

### Same but without length Norm

Chimp.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Chimp$Counts[1:NumRowsToAnalyze,], rep(1, length(CountTables$Chimp$GeneLengths[1:NumRowsToAnalyze])))
Human.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Human$Counts[1:NumRowsToAnalyze,], rep(1, length(CountTables$Chimp$GeneLengths[1:NumRowsToAnalyze])))

cbind(Chimp.NB.fit.parameters,
      apply(CountTables$Chimp$log2RPKM, 1, mean)[1:NumRowsToAnalyze],
      Human.NB.fit.parameters,
      apply(CountTables$Human$log2RPKM, 1, mean)[1:NumRowsToAnalyze]) %>%
  setNames(
    c("Chimp.Mean.Expression", "Chimp.Overdispersion", "Chimp.theta.se", "Chimp.Mean.Log2RPKM", "Human.Mean.Expression", "Human.Overdispersion", "Human.theta.se","Human.Mean.Log2RPKM")
  ) %>%
  rownames_to_column(var="gene") %>%
  mutate(
    Chimp.Residual = GetLoessResidual(Chimp.Mean.Expression, log(Chimp.Overdispersion)),
    Human.Residual = GetLoessResidual(Human.Mean.Expression, log(Human.Overdispersion))

  ) %>%
  mutate(
    Chimp.Dispersion = Chimp.Residual + median(log(Chimp.Overdispersion), na.rm = T),
    Human.Dispersion = Human.Residual + median(log(Human.Overdispersion), na.rm = T)
  ) %>%
  dplyr::select(gene, everything()) %>%
  write.table("../../output/OverdispersionEstimatesFromChimp.NoLengthNorm.txt", sep='\t', quote=F, row.names = F)


### Estimate overdispesion but after dropping virus challanged chimps
CountTables <- GetCountTables(CountTableChimpFile,
                              CountTableHumanFile,
                              0, GeneListForOverdispersionCalculation, ChimpSampleDrop=c(ChimpSamplesToDrop, VirusChallengedChimps), HumanSampleDrop = HumanSamplesToDrop)


Chimp.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Chimp$Counts[1:NumRowsToAnalyze,], CountTables$Chimp$GeneLengths[1:NumRowsToAnalyze])
Human.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Human$Counts[1:NumRowsToAnalyze,], CountTables$Human$GeneLengths[1:NumRowsToAnalyze])

cbind(Chimp.NB.fit.parameters,
      apply(CountTables$Chimp$log2RPKM, 1, mean)[1:NumRowsToAnalyze],
      Human.NB.fit.parameters,
      apply(CountTables$Human$log2RPKM, 1, mean)[1:NumRowsToAnalyze]) %>%
  setNames(
    c("Chimp.Mean.Expression", "Chimp.Overdispersion", "Chimp.theta.se", "Chimp.Mean.Log2RPKM", "Human.Mean.Expression", "Human.Overdispersion", "Human.theta.se","Human.Mean.Log2RPKM")
  ) %>%
  rownames_to_column(var="gene") %>%
  mutate(
    Chimp.Residual = GetLoessResidual(Chimp.Mean.Expression, log(Chimp.Overdispersion)),
    Human.Residual = GetLoessResidual(Human.Mean.Expression, log(Human.Overdispersion))

  ) %>%
  mutate(
    Chimp.Dispersion = Chimp.Residual + median(log(Chimp.Overdispersion), na.rm = T),
    Human.Dispersion = Human.Residual + median(log(Human.Overdispersion), na.rm = T)
  ) %>%
  dplyr::select(gene, everything()) %>%
  write.table("../../output/OverdispersionEstimatesFromChimp_NoVirusChallangedIndividuals.txt", sep='\t', quote=F, row.names = F)

### Same but without length Norm

Chimp.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Chimp$Counts[1:NumRowsToAnalyze,], rep(1, length(CountTables$Chimp$GeneLengths[1:NumRowsToAnalyze])))
Human.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Human$Counts[1:NumRowsToAnalyze,], rep(1, length(CountTables$Chimp$GeneLengths[1:NumRowsToAnalyze])))

cbind(Chimp.NB.fit.parameters,
      apply(CountTables$Chimp$log2RPKM, 1, mean)[1:NumRowsToAnalyze],
      Human.NB.fit.parameters,
      apply(CountTables$Human$log2RPKM, 1, mean)[1:NumRowsToAnalyze]) %>%
  setNames(
    c("Chimp.Mean.Expression", "Chimp.Overdispersion", "Chimp.theta.se", "Chimp.Mean.Log2RPKM", "Human.Mean.Expression", "Human.Overdispersion", "Human.theta.se","Human.Mean.Log2RPKM")
  ) %>%
  rownames_to_column(var="gene") %>%
  mutate(
    Chimp.Residual = GetLoessResidual(Chimp.Mean.Expression, log(Chimp.Overdispersion)),
    Human.Residual = GetLoessResidual(Human.Mean.Expression, log(Human.Overdispersion))

  ) %>%
  mutate(
    Chimp.Dispersion = Chimp.Residual + median(log(Chimp.Overdispersion), na.rm = T),
    Human.Dispersion = Human.Residual + median(log(Human.Overdispersion), na.rm = T)
  ) %>%
  dplyr::select(gene, everything()) %>%
  write.table("../../output/OverdispersionEstimatesFromChimp_NoVirusChallangedIndividuals.NoLengthNorm.txt", sep='\t', quote=F, row.names = F)

