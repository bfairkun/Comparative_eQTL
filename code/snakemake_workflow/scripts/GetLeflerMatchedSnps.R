library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
plink_ld_out <- args[1]
snps_location_out <- args[2]

# plink_ld_out<-"MiscOutput/SpeciesSharedSNPs.ld"
# setwd("/project/gilad/bjf79/projects/Comparative_eQTL/code/snakemake_workflow/")

LD <- read.table(plink_ld_out, header=TRUE)
# head(LD)

# LD %>% distinct(SNP_A, .keep_all = T) %>% dim()

LD %>% sample_frac() %>%
  filter(abs(MAF_A-MAF_B)<=0.05) %>%
  filter(R2<=0.2) %>%
  distinct(SNP_B, .keep_all = T) %>%
  distinct(SNP_A, .keep_all = T) %>%
  select(SNP_B,CHR_B, BP_B) %>%
  write.table(snps_location_out, row.names=F, col.names = c("snp", "chrom", "pos"), sep='\t', quote=F)
