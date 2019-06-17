library(tidyverse)

SNP_file_name <- args[1]
MatrixEqtlResults <- args[2]

MatrixEqtlResults <- "output/MatrixEQTL_sig_eqtls.txt"
SNP_file_name <- "output/MatrixEQTL_sig_genotypes.raw"
PhenotypeTable <- "output/log10TPM.StandardizedAndNormalized.txt"
FamFile <- "output/ForAssociationTesting.temp.fam"
DirectoryOut <- "code/snakemake_workflow/scratch/boxplots/"

eQTLs <- read.table(MatrixEqtlResults, header=T)
SampleList <- read.table(FamFile)$V2

Genotypes <- read.table(SNP_file_name, header=T, check.names = F, stringsAsFactors = F)
colnames(Genotypes) <- sub("_.*", "", colnames(Genotypes))

Phenotypes <- read.table(PhenotypeTable, header=T, check.names=FALSE, row.names = 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "FID") %>%
  filter(FID %in% SampleList)

MergedData <- left_join(Genotypes, Phenotypes, by=c("IID" = "FID"))  %>%
  as.data.frame()

MyBoxplot <- function(DataFrame, Labels.name, SNP.name, Gene.name){
  data.frame(Genotype = DataFrame[[SNP.name]],
             Phenotype = DataFrame[[Gene.name]],
             FID=DataFrame[[Labels.name]]) %>%
    ggplot(aes(x=factor(Genotype), y=Phenotype, label=FID)) +
    geom_boxplot(outlier.shape = NA) +
    geom_text(position=position_jitter(width=0.25), alpha=1, size=2) +
    scale_y_continuous(name=paste(Gene.name, "\nnormalized expression")) +
    xlab(paste(SNP.name, "Genotype")) +
    theme_bw()
}

BestSnps <- eQTLs %>%
  group_by(gene) %>%
  dplyr::slice(which.min(qvalue))

for (i in 1:nrow(BestSnps)){
  try({
    A<-MyBoxplot(MergedData, "IID", as.character(unlist(BestSnps[i,"snps"])), as.character(unlist(BestSnps[i,"gene"])))
    ggsave(filename = paste0(DirectoryOut, as.character(unlist(BestSnps[i,"gene"])), ".png"), plot=A)
  })
}
