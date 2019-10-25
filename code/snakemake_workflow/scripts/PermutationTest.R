library(tidyverse)
library(data.table)
library(qvalue)
library(stats)


args = commandArgs(trailingOnly=TRUE)
PermutationsTableFile <- args[1]
ActualTableFile <- args[2]
TableOutFile <- args[3]

## Files for testing
# setwd("~/CurrentProjects/sQTL_pipeline/")
# PermutationsTableFile <- "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/PermutationsCombined.txt"
# ActualTableFile <- "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/BestPvalueNonPermuted.txt"
# TableOutFile <- "~/temp/test.txt"

# read table of min p values from permutations
permutationTable <- t(fread(PermutationsTableFile, sep='\t', header=T))

# read table of min p values from real data
ActualTable <- fread(ActualTableFile, sep='\t', header=T) %>%
  as.data.frame() %>% tibble::column_to_rownames("Gene")

# Merge the two tables, putting the real data as the first column.
# If you are interested in gene level or cluster-level tests, this where you
# should group and find minimum Pval. Here I did cluster level P-values.
MergedTable <- merge(ActualTable, permutationTable, by="row.names") %>%
  mutate(clusterName=gsub("^(.+?):.+?:.+?:(clu_\\d+).+$", "\\1.\\2", Row.names, perl=T)) %>%
  dplyr::group_by(clusterName) %>% dplyr::summarise_all(dplyr::funs(min)) %>%
  dplyr::ungroup() %>% dplyr::select(-Row.names) %>% as.data.frame() %>%
  tibble::column_to_rownames("clusterName")

PermutationTest <- function(PvalVector){
  return(ecdf(PvalVector)(PvalVector[1]))
}


OutTable <- data.frame(
  BestActualPval=MergedTable$MinP,
  PermutationTestPval=apply(MergedTable, 1, PermutationTest))
OutTable$PermutationTestQvalue <- qvalue(OutTable$PermutationTestPval)$qvalue

OutTable %>%
  tibble::rownames_to_column(var = "GroupedTrait") %>%
  write.table(file=TableOutFile, quote=F, sep='\t')
