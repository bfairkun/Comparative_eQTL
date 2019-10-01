
# functions to source in later analysis to make near publication ready figures

# function to read in four files and return list of dataframes which often used later:
# - Human tsv file of gene (sysID), pvalue, beta
# - ...same for chimp
# - SysToID Biomart file that contains dN/dS and other info

library(tidyverse)
library(MASS)
library(edgeR)

# return matrix of column x repeated n times
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


# Read in relevant files of Pvalues for all genes tested and the SysToID file that I downloaded from Biomart
function(Chimp.tsv, Human.tsv, SysToID.tsv){

}


# Return plots for dN/dS, etc. from dataframe
function(Chimp.df, Human.df, SysToID.df){

}

# Return intercept of NB fit (estimate of mean expression our context), given Y (a vector of counts for a single gene), and libsizes (total counts for libraries)
GetInterceptOfNB.fit <- function(Y, libsizes, genesize=1){
  tryCatch(
    expr = {
      fit <- MASS::glm.nb(Y ~ offset(log(libsizes))-offset(log(GeneSizeVector))  + 1)
      intercept <- fit$coefficients
    },
    error = function(e){intercept <- NA},
    finally = {
      return(intercept)
    }
  )
}

# Return 1/theta of NB fit (estimate of overdispersion in our context), given Y (a vector of counts for a single gene), and libsizes (total counts for libraries), and optionally, vector of gene sizes
GetOverdispersionOfUnderlyingGamma <- function(Y, libsizes, genesize=1){
  tryCatch(
    expr = {
      GeneSizeVector <- rep(genesize, length(libsizes))
      fit <- MASS::glm.nb(Y ~ offset(log(libsizes))-offset(log(GeneSizeVector))  + 1)
      # var <- fit$coefficients**2/fit$theta
      Overdispersion <- 1/fit$theta

    },
    error = function(e){Overdispersion <- NA},
    finally = {
      return(Overdispersion)
    }
  )
}

# Return list of 1/theta and intercept of NB.fit call from Abhishek's snippet
Get.NB.Fit.Parameters <- function(Y, libsizes, genesize=1){
  tryCatch(
    {
      GeneSizeVector <- rep(genesize, length(libsizes))
      fit <- MASS::glm.nb(Y ~ offset(log(libsizes))-offset(log(GeneSizeVector))  + 1)
      # overdispersion <- 1/fit$theta
      # intercept <- fit$coefficients
      return(list(Overdispersion=1/fit$theta, Intercept=fit$coefficients))

    },
    error = function(e){
      return(list(Overdispersion=NA, Intercept=NA))
      }
  )
}

# Return dataframe of overdispersion estimates and mu, given a count table
GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable <- function(MyCountTable, genesizes){
  out <- matrix(data=NA, nrow=nrow(MyCountTable), ncol=2)
  LibSize <- colSums(MyCountTable)
  for(i in 1:nrow(MyCountTable)) {
    NB.fit <- Get.NB.Fit.Parameters(t(MyCountTable[i,]), LibSize, genesizes[i])
    out[i,1] <- NB.fit$Intercept
    out[i,2] <- NB.fit$Overdispersion
    rownames(out) <- rownames(MyCountTable)
    colnames(out) <- c("mu", "overdispersion")
  }
  return(as.data.frame(out))
}



# Return a named vector of standard dev independent of log(mean) trend from a count table
function(Count.table.df){

}

# Return a named vector of overdispersion independent of log(mean) trend from a count table
function(Count.table.df){

}

# Return a plot and GSEA analysis
function(NamedAndOrderedVectorOfEnsemblGenes){

}

GetCountTables <-function(ChimpCountTableFile, HumanCountTableFile, SubsampleSize, GenesToKeep, ChimpSampleDrop=NULL, HumanSampleDrop=NULL)
  #if SubsampleSize parameter == 0, use full table, otherwise, subsample from it
{
  FullChimpData <- read.table(gzfile(ChimpCountTableFile), header=T, check.names=FALSE, skip=1)
  FullHumanData <- read.table(gzfile(HumanCountTableFile), header=T, check.names=FALSE, skip=1)

  if (!is.null(ChimpSampleDrop)){
    FullChimpData <- FullChimpData %>% dplyr::select(-ChimpSampleDrop)
  }
  if (!is.null(HumanSampleDrop)){
    FullHumanData <- FullHumanData %>% dplyr::select(-HumanSampleDrop)
  }

  if (SubsampleSize==0){
    CountTableChimp <- FullChimpData
    colnames(CountTableChimp) <- paste0("C.", colnames(CountTableChimp))
    CountTableHuman <- FullHumanData
    colnames(CountTableHuman) <- paste0("H.", colnames(CountTableHuman))

  } else {
    CountTableChimp <- FullChimpData %>% dplyr::select(c(1:6, sample(7:length(FullChimpData), SubsampleSize)))
    colnames(CountTableChimp) <- paste0("C.", colnames(CountTableChimp))

    CountTableHuman <- FullHumanData %>% dplyr::select(c(1:6, sample(7:length(FullHumanData), SubsampleSize)))
    colnames(CountTableHuman) <- paste0("H.", colnames(CountTableHuman))
  }

  CombinedTable <- inner_join(CountTableChimp[,c(1,7:length(CountTableChimp))], CountTableHuman[,c(1,7:length(CountTableHuman))], by=c("C.Geneid"="H.Geneid")) %>%
    column_to_rownames("C.Geneid") %>% as.matrix()

  SpeciesFactor <- colnames(CombinedTable) %>% substr(1,1) %>% factor() %>% unclass() %>% as.character()

  d0 <- DGEList(CombinedTable)
  d0 <- calcNormFactors(d0)



  d <- d0[GenesToKeep,]
  mm <- model.matrix(~0 + SpeciesFactor)
  cpm.table <- cpm(d, log=T)

  GeneLengths <- left_join((d$counts %>% as.data.frame() %>% rownames_to_column() %>% dplyr::select(rowname)),
    CountTableChimp[,c("C.Geneid", "C.Length")], by=c("rowname"="C.Geneid")) %>%
    left_join(CountTableHuman[,c("H.Geneid", "H.Length")], by=c("rowname"="H.Geneid"))
  GeneLengthMatrix <- cbind(
    rep.col(log2(GeneLengths$C.Length/1000), length(CountTableChimp)-6),
    rep.col(log2(GeneLengths$H.Length/1000), length(CountTableHuman)-6))
  rownames(GeneLengthMatrix) <- GeneLengths$rowname
  log2RPKM.table <- cpm.table - GeneLengthMatrix[rownames(cpm.table),]

  # separate the count tables into human and chimp tables
  ChimpLog2RPKM.table <- log2RPKM.table %>% as.data.frame() %>%
    dplyr::select(contains("C."))
  HumanLog2RPKM.table <- log2RPKM.table %>% as.data.frame() %>%
    dplyr::select(contains("H."))
  ChimpCount.table <- d$counts %>% as.data.frame() %>%
    dplyr::select(contains("C."))
  HumanCount.table <- d$counts %>% as.data.frame() %>%
    dplyr::select(contains("H."))

  ToReturn <- list(
    Chimp=list(Counts=ChimpCount.table, log2RPKM=ChimpLog2RPKM.table, GeneLengths=GeneLengths$C.Length),
    Human=list(Counts=HumanCount.table, log2RPKM=HumanLog2RPKM.table, GeneLengths=GeneLengths$H.Length)
  )


  return(ToReturn)
}


### tests


# setwd("analysis/")
# HumanSamplesToDrop <- c(c("SRR1507229","SRR603918", "SRR1478149", "SRR598509", "SRR613186"), c("SRR1489693", "SRR598148", "59167", "SRR1478900", "SRR1474730", "61317"))
# ChimpSamplesToDrop <- c("Little_R")
# CountTableChimpFile <- '../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz'
# CountTableHumanFile <- '../output/PowerAnalysisFullCountTable.Human.subread.txt.gz'
# eQTLs <- read.table(gzfile("../data/PastAnalysesDataToKeep/20190521_eQTLs_250kB_10MAF.txt.gz"), header=T)
# # List of chimp tested genes
# ChimpTestedGenes <- rownames(read.table('../output/ExpressionMatrix.un-normalized.txt.gz', header=T, check.names=FALSE, row.names = 1))
#
# ChimpToHumanGeneMap <- read.table("../data/Biomart_export.Hsap.Ptro.orthologs.txt.gz", header=T, sep='\t', stringsAsFactors = F)
#
# # Of this ortholog list, how many genes are one2one
# OneToOneMap <- ChimpToHumanGeneMap %>%
#   filter(Chimpanzee.homology.type=="ortholog_one2one")
#
# # Read gtex heart egene list
# # Only consider those that were tested in both species and are one2one orthologs
# GtexHeartEgenes <- read.table("../data/Heart_Left_Ventricle.v7.egenes.txt.gz", header=T, sep='\t', stringsAsFactors = F) %>%
#   mutate(gene_id_stable = gsub(".\\d+$","",gene_id)) %>%
#   filter(gene_id_stable %in% OneToOneMap$Gene.stable.ID) %>%
#   mutate(chimp_id = plyr::mapvalues(gene_id_stable, OneToOneMap$Gene.stable.ID,  OneToOneMap$Chimpanzee.gene.stable.ID, warn_missing = F)) %>%
#   filter(chimp_id %in% ChimpTestedGenes)
#
#
# EgenesTested <- gsub("\\..+", "", GtexHeartEgenes$gene_id, perl=T)
# length(EgenesTested)
# GenesInDESet <- read.table(gzfile(CountTableChimpFile), header=T, check.names=FALSE, skip=1)$Geneid
# length(GenesInDESet)
#
# GeneList <- intersect(as.character(GenesInDESet),EgenesTested)
# length(GeneList)
#
#
# CountTables <- GetCountTables(CountTableChimpFile,
#                                CountTableHumanFile,
#                                0, GeneList, ChimpSampleDrop=ChimpSamplesToDrop, HumanSampleDrop = HumanSamplesToDrop)
#
# NumRowsToAnalyze=5000
# A<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Chimp$Counts[1:NumRowsToAnalyze,], CountTables$Chimp$GeneLengths[1:NumRowsToAnalyze])
# B<-GetParameterEstimatesOfUnderlyingGamma_lengthAdjusted_FromTable(CountTables$Human$Counts[1:NumRowsToAnalyze,], CountTables$Human$GeneLengths[1:NumRowsToAnalyze])
# qplot(A$overdispersion,B$overdispersion) +
#   scale_x_continuous(trans='log10', limits=c(1E-2,10)) +
#   scale_y_continuous(trans='log10', limits=c(1E-2,10)) +
#   theme_bw()
#
# cor(A$mu,log(A$overdispersion), use="complete.obs")
# qplot(A$mu, A$overdispersion, alpha=0.05) +
#   scale_y_continuous(trans='log10', limits=c(1E-2,10)) +
#   theme_bw()
#
# cor(A$mu,log(A$overdispersion), use="complete.obs")
# qplot(apply(CountTables$Chimp$log2RPKM[1:NumRowsToAnalyze,],1,mean), sqrt(apply(CountTables$Chimp$log2RPKM[1:NumRowsToAnalyze,],1,var)), alpha=0.05) +
#   theme_bw()

###

