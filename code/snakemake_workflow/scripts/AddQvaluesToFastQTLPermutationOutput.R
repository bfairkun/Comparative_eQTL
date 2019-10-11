library(tidyverse)
library(qvalue)

## FastQTL outputs 10 fields:
# ID of the tested molecular phenotype (in this particular case, the gene ID)
# Number of variants tested in cis for this phenotype
# MLE of the shape1 parameter of the Beta distribution
# MLE of the shape2 parameter of the Beta distribution
# Dummy [To be described later]
# ID of the best variant found for this molecular phenotypes (i.e. with the smallest p-value)
# Distance between the molecular phenotype - variant pair
# The nominal p-value of association that quantifies how significant from 0, the regression coefficient is
# The slope associated with the nominal p-value of association [only in version > v2-184]
# A first permutation p-value directly obtained from the permutations with the direct method. This is basically a corrected version of the nominal p-value that accounts for the fact that multiple variants are tested per molecular phenotype.
# A second permutation p-value obtained via beta approximation. We advice to use this one in any downstream analysis.
## This script will add storey's qvalue as a 11th field

args = commandArgs(trailingOnly=TRUE)
FastQTLIn <- args[1]
Output <- args[2]

# setwd("~/CurrentProjects/Comparative_eQTL/code/snakemake_workflow/")
# FastQTLIn <- "GTEX_renalysis/SubsamplesVaryingSize/output/SampleSize_200.txt.gz"
# Output <- "scratch/TestFastQTL.AddQvalues.txt"

d <- read.table(FastQTLIn, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
d$st = qvalue(d$bpval)$qvalues
write.table(d, Output, quote=F, sep='\t', row.names = F)


