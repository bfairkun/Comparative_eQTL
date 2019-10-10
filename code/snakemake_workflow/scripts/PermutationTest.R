library(tidyverse)
library(data.table)
library(qvalue)
library(stats)

setwd("/project2/gilad/bjf79_project1/projects/Comparative_eQTL/code/snakemake_workflow/")

permutationTable <- t(fread("eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/PermutationsCombined.txt", sep='\t', header=T))
ActualTable <- fread("eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt", sep='\t', header=T)
eigenMTTable <- read.table("eQTL_mapping/MatrixEQTL/eigenMT.corrected.eGenes.txt.gz", header=T, sep='\t') %>%
  tibble::column_to_rownames("gene")


BestActualPVal <- ActualTable %>%
  dplyr::group_by(gene) %>%
  dplyr::slice(which.min(pvalue)) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, pvalue) %>%
  tibble::column_to_rownames("gene") %>% as.matrix()

MergedTable <- merge(BestActualPVal, permutationTable, by="row.names") %>%
  tibble::column_to_rownames("Row.names")

PermutationTest <- function(PvalVector){
  return(ecdf(PvalVector)(PvalVector[1]))
}

PermutationTestPvals <- apply(MergedTable, 1, PermutationTest)

B<-merge(as.data.frame(BestActualPVal), as.data.frame(PermutationTestPvals), by="row.names") %>%
rownames(B) <- B$Row.names
C<-merge(B, eigenMTTable, by="row.names")

plot(-log10(C$PermutationTestPvals), -log10(C$FDR))

