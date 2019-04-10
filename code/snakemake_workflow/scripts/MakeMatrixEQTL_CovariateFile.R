library(plyr)
library(tidyverse)
library(readxl)
library(data.table)

GenotypePCs <- read.table("output/PopulationStructure/pca.eigenvec", header=T) %>%
  select(IID, PC1, PC2, PC3) %>%
  dplyr::rename(Individual.ID=IID, GenotypePC1=PC1, GenotypePC2=PC2, GenotypePC3=PC3)

OtherMetadata <- as.data.frame(read_excel("data/Metadata.xlsx")) %>%
  select(Individual.ID, SX) %>%
  transmute(Individual.ID=Individual.ID, SX = plyr::mapvalues(SX, c("M", "F"), c(0,1)))



fam_file_name = "code/snakemake_workflow/eQTL_mapping/plink/ForAssociationTesting.fam"
expression<-fread(fam_file_name, header=F)

PCResults <- expression %>%
  select(-c(V1,V2,V3,V4,V5)) %>%
  prcomp(center=T, scale. = T)

as.data.frame(PCResults$x, row.names = expression$V1) %>%
  rownames_to_column(var = "Individual.ID") %>%
  select(Individual.ID, PC1:PC20) %>%
  left_join(OtherMetadata, by="Individual.ID") %>%
  left_join(GenotypePCs, by= "Individual.ID") %>%
  select(Individual.ID, SX, GenotypePC1, GenotypePC2, GenotypePC3, PC1, PC2, PC3, PC4, PC5, PC6, PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15) %>%
  t() %>%
  write.table(file="output/Covariates.15.txt", sep='\t', quote = F, row.names = T, col.names = F)
