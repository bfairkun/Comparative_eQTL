library(data.table)
library(tidyverse)
GTEXCounts <- fread(paste0("gunzip -c ", "/project/gilad/bjf79/projects/Comparative_eQTL/code/snakemake_workflow/scratch/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz"), sep="\t", header=TRUE, data.table=FALSE)
HeartSamples <- read.table("/project/gilad/bjf79/projects/Comparative_eQTL/code/snakemake_workflow/scratch/GtexHeartSampleList.txt")

GTEXCounts %>% select(HeartSamples$V1) %>%
  write.table(file="/project/gilad/bjf79/projects/Comparative_eQTL/code/snakemake_workflow/scratch/GTEx_HeartCountTable.txt", sep="\t", row.names = F, quote=F)
