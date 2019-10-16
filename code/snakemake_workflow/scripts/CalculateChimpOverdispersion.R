# setwd("../")

source("../CustomFunctions.R")
library(readxl)

HumanSamplesToDrop <- c(c("SRR1507229","SRR603918", "SRR1478149", "SRR598509", "SRR613186"), c("SRR1489693", "SRR598148", "59167", "SRR1478900", "SRR1474730", "61317"))
ChimpSamplesToDrop <- c("Little_R")
OtherMetadata <- as.data.frame(read_excel("../../data/Metadata.xlsx"))
VirusChallengedChimps <- OtherMetadata %>% filter(grepl("V+",Viral.status)) %>% pull(Individual.ID)

CountTableChimpFile <- '../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz'
CountTableHumanFile <- '../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz'

EgenesTested <- TsvToCombinedEgenes(Chimp.tsv = "../../output/ChimpEgenes.eigenMT.txt.gz", Human.tsv = "../../data/GTEX_v8_eGenes/Heart_Left_Ventricle.v8.egenes.txt.gz", SysToID.tsv = "../../data/Biomart_export.Hsap.Ptro.orthologs.txt.gz", HumanTsvType = "GTEx")
GeneListTestedForQTLs <- EgenesTested$H.gene
GenesInDESet <- read.table(gzfile(CountTableChimpFile), header=T, check.names=FALSE, skip=1)$Geneid
GeneListForOverdispersionCalculation <- intersect(as.character(GenesInDESet),GeneListTestedForQTLs)



### estimate overdispersion
CountTables <- GetCountTables(CountTableChimpFile,
                              CountTableHumanFile,
                              0, GeneListForOverdispersionCalculation, ChimpSampleDrop=ChimpSamplesToDrop, HumanSampleDrop = HumanSamplesToDrop)

NumRowsToAnalyze=length(GeneListForOverdispersionCalculation)
# NumRowsToAnalyze=10


Chimp.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Chimp$Counts[1:NumRowsToAnalyze,], CountTables$Chimp$GeneLengths[1:NumRowsToAnalyze])
Human.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Human$Counts[1:NumRowsToAnalyze,], CountTables$Human$GeneLengths[1:NumRowsToAnalyze])

ToPlot <- cbind(Chimp.NB.fit.parameters,
                apply(CountTables$Chimp$log2RPKM, 1, mean),
                Human.NB.fit.parameters,
                apply(CountTables$Human$log2RPKM, 1, mean))
colnames(ToPlot) <- c("Chimp.Mean.Expression", "Chimp.Overdispersion", "Chimp.theta.se", "Chimp.Mean.Log2RPKM", "Human.Mean.Expression", "Human.Overdispersion", "Human.theta.se","Human.Mean.Log2RPKM")

GetLoessResidual <- function(x, y){
  loess.fit <- loess(y ~ x, degree=1)
  #for our purposes, the fit with the default degree=2 looks like unnatural fit, especially near the edges
  loess.residual <- y - predict(loess.fit, x)
  return(loess.residual)
}

ToPlot$Chimp.Residual <- exp(GetLoessResidual(ToPlot$Chimp.Mean.Expression, log(ToPlot$Chimp.Overdispersion)))
ToPlot$Human.Residual <- exp(GetLoessResidual(ToPlot$Human.Mean.Expression, log(ToPlot$Human.Overdispersion)))

ToPlot %>%
  rownames_to_column(var="gene") %>%
  dplyr::select(gene, everything()) %>%
write.table("../../output/OverdispersionEstimatesFromChimp.txt", sep='\t', quote=F, row.names = F)


### Estimate overdispesion but after dropping virus challanged chimps
CountTables <- GetCountTables(CountTableChimpFile,
                              CountTableHumanFile,
                              0, GeneListForOverdispersionCalculation, ChimpSampleDrop=c(ChimpSamplesToDrop, VirusChallengedChimps), HumanSampleDrop = HumanSamplesToDrop)

NumRowsToAnalyze=length(GeneListForOverdispersionCalculation)
# NumRowsToAnalyze=10


Chimp.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Chimp$Counts[1:NumRowsToAnalyze,], CountTables$Chimp$GeneLengths[1:NumRowsToAnalyze])
Human.NB.fit.parameters<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Human$Counts[1:NumRowsToAnalyze,], CountTables$Human$GeneLengths[1:NumRowsToAnalyze])

ToPlot <- cbind(Chimp.NB.fit.parameters,
                apply(CountTables$Chimp$log2RPKM, 1, mean),
                Human.NB.fit.parameters,
                apply(CountTables$Human$log2RPKM, 1, mean))
colnames(ToPlot) <- c("Chimp.Mean.Expression", "Chimp.Overdispersion", "Chimp.theta.se", "Chimp.Mean.Log2RPKM", "Human.Mean.Expression", "Human.Overdispersion", "Human.theta.se","Human.Mean.Log2RPKM")

GetLoessResidual <- function(x, y){
  loess.fit <- loess(y ~ x, degree=1)
  #for our purposes, the fit with the default degree=2 looks like unnatural fit, especially near the edges
  loess.residual <- y - predict(loess.fit, x)
  return(loess.residual)
}

ToPlot$Chimp.Residual <- exp(GetLoessResidual(ToPlot$Chimp.Mean.Expression, log(ToPlot$Chimp.Overdispersion)))
ToPlot$Human.Residual <- exp(GetLoessResidual(ToPlot$Human.Mean.Expression, log(ToPlot$Human.Overdispersion)))

ToPlot %>%
  rownames_to_column(var="gene") %>%
  dplyr::select(gene, everything()) %>%
  write.table("../../output/OverdispersionEstimatesFromChimp_NoVirusChallangedIndividuals.txt", sep='\t', quote=F, row.names = F)
