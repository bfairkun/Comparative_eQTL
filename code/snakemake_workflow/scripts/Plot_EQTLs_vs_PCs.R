library(tidyverse)
library(gridExtra)
library(reshape2)

args <- commandArgs(trailingOnly = T)
FileIn <- args[1]
PlotOut <- args[2]

results <- read.table(FileIn, sep='\t', header=F, col.names = c("SNP", "gene", "beta","t", "p", "q", "Filename"))
head(results)
## eGene count
eGenes_perPC <- results %>%
  filter(q<0.2) %>%
  distinct(gene, Filename, .keep_all = TRUE) %>%
  group_by(Filename) %>%
  tally() %>%
  mutate(PCcount = gsub("^(.+?)_and_(.+?)RNASeqPCs.covariates.txt","\\2",Filename, perl=TRUE)) %>%
  select(-Filename, FDR_20=n)

eGenes_perPC$FDR_10 <- results %>%
  filter(q<0.1) %>%
  distinct(gene, Filename, .keep_all = TRUE) %>%
  group_by(Filename) %>%
  tally() %>%
  mutate(PCcount = gsub("^(.+?)_and_(.+?)RNASeqPCs.covariates.txt","\\2",Filename, perl=TRUE)) %>%
  pull(n)

eGenes_perPC$FDR_5 <- results %>%
  filter(q<0.05) %>%
  distinct(gene, Filename, .keep_all = TRUE) %>%
  group_by(Filename) %>%
  tally() %>%
  mutate(PCcount = gsub("^(.+?)_and_(.+?)RNASeqPCs.covariates.txt","\\2",Filename, perl=TRUE)) %>%
  pull(n)

P1 <- eGenes_perPC %>%
  melt(id.vars="PCcount") %>%
  transform(PCcount=as.numeric(PCcount)) %>%
  mutate(FDR = plyr::mapvalues(variable, from=c("FDR_10", "FDR_20", "FDR_5"), to=c("0.1", "0.2", "0.05"))) %>%
  ggplot(aes(x=PCcount, y=value, color=FDR)) +
  geom_point() +
  geom_line() +
  xlab("Number PCs") +
  ylab("Number eGenes") +
  theme_bw()

## eQTL count
eGenes_perPC <- results %>%
  filter(q<0.2) %>%
  group_by(Filename) %>%
  tally() %>%
  mutate(PCcount = gsub("^(.+?)_and_(.+?)RNASeqPCs.covariates.txt","\\2",Filename, perl=TRUE)) %>%
  select(-Filename, FDR_20=n)

eGenes_perPC$FDR_10 <- results %>%
  filter(q<0.1) %>%
  group_by(Filename) %>%
  tally() %>%
  mutate(PCcount = gsub("^(.+?)_and_(.+?)RNASeqPCs.covariates.txt","\\2",Filename, perl=TRUE)) %>%
  pull(n)

eGenes_perPC$FDR_5 <- results %>%
  filter(q<0.05) %>%
  group_by(Filename) %>%
  tally() %>%
  mutate(PCcount = gsub("^(.+?)_and_(.+?)RNASeqPCs.covariates.txt","\\2",Filename, perl=TRUE)) %>%
  pull(n)

P2 <- eGenes_perPC %>%
  melt(id.vars="PCcount") %>%
  transform(PCcount=as.numeric(PCcount)) %>%
  mutate(FDR = plyr::mapvalues(variable, from=c("FDR_10", "FDR_20", "FDR_5"), to=c("0.1", "0.2", "0.05"))) %>%
  ggplot(aes(x=PCcount, y=value, color=FDR)) +
  geom_point() +
  geom_line() +
  xlab("Number PCs") +
  ylab("Number eGene-SNP pairs") +
  theme_bw()

g<-arrangeGrob(P1, P2, nrow=2)
ggsave(PlotOut, g)
