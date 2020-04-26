library(tidyverse)
library(qvalue)

args = commandArgs(trailingOnly=TRUE)
DispersionF <- args[1]
DispersionSE_F <- args[2]
DispersionP_F <- args[3]
outputFn <- args[4]

Dispersion <- read.delim(DispersionF)
Dispersion.se <- read.delim(DispersionSE_F)
Dispersion.P <- read.delim(DispersionP_F)

GeneNames <- read.delim('../../data/HumanGeneIdToHGNC_Symbol.Biomart.txt.gz')


TableS3 <- Dispersion %>%
  dplyr::select(gene,
                Chimp.Mean.Expression,
                Human.Mean.Expression,
                Chimp.Overdispersion,
                Human.Overdispersion,
                Chimp.Mean.Adjusted.Dispersion=Chimp.Residual,
                Human.Mean.Adjusted.Dispersion=Human.Residual) %>%
  left_join(GeneNames, by=c("gene"="Gene.stable.ID")) %>%
  left_join(Dispersion.se, by="gene") %>%
  left_join(Dispersion.P, by="gene")
hist(TableS3$P)
qobj <- qvalue(p = TableS3$P)
table(qobj$qvalues < 0.1)
TableS3 %>%
  mutate(q.value=qobj$qvalues) %>%
  dplyr::select(gene, HGNC.symbol, everything()) %>%
  dplyr::select(-HGNC.ID) %>%
  write_delim(outputFn, delim='\t')
