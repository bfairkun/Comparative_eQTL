library(tidyverse)
library("edgeR")
source("../CustomFunctions.R")

CountTableChimpFile <- '../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz'
CountTableHumanFile <- '../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz'
OutputDE <- '../../output/Final/TableS1.tab'
GeneIDs <- '../../data/HumanGeneIdToHGNC_Symbol.Biomart.txt.gz'

HumanSamplesToDrop <- c("SRR613186",  "SRR598509",  "SRR1478149", "SRR603918",  "SRR1507229", "SRR1478900", "SRR1477015", "SRR601986",  "SRR614996",  "SRR1474730") 
ChimpSamplesToDrop <- c()

GeneID.df <- read.delim(GeneIDs) %>%
  dplyr::select(-HGNC.ID)

CountTableChimp <- read.table(gzfile(CountTableChimpFile), header=T, check.names=FALSE, skip=1)
colnames(CountTableChimp) <- paste0("C.", colnames(CountTableChimp))

CountTableHuman <- read.table(gzfile(CountTableHumanFile), header=T, check.names=FALSE, skip=1)
colnames(CountTableHuman) <- paste0("H.", colnames(CountTableHuman))

CombinedTable <- inner_join(CountTableChimp[,c(1,7:length(CountTableChimp))], CountTableHuman[,c(1,7:length(CountTableHuman))], by=c("C.Geneid"="H.Geneid")) %>%
  column_to_rownames("C.Geneid") %>%
  dplyr::select(-c(paste0("H.", HumanSamplesToDrop))) %>%
  as.matrix()

cpm <- cpm(CombinedTable, log=TRUE, prior.count=0.5)
d0 <- DGEList(CombinedTable)

#Calculate normalization factors
d0 <- calcNormFactors(d0)
cutoff <- 6
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)

# DE.Subsampled <-function(ChimpCountTableFile, HumanCountTableFile, SubsampleSize, FC.NullInterval, drop, ChimpSampleDrop=NULL, HumanSampleDrop=NULL)
FC.NullInterval <- log2(1.0)

True.efit <- DE.Subsampled(CountTableChimpFile, CountTableHumanFile, SubsampleSize = 39, FC.NullInterval = FC.NullInterval, drop=drop, ChimpSampleDrop = NULL, HumanSampleDrop = HumanSamplesToDrop)

TrueResponse <- decideTests(True.efit, p.value=0.05)
TableOut <- toptable(True.efit, n=Inf, sort.by="none") %>%
  rownames_to_column("Ensembl_geneID") %>%
  left_join(GeneID.df, by=c("Ensembl_geneID"="Gene.stable.ID")) %>%
  dplyr::select(Ensembl_geneID,
                HGNC.symbol,
                logFC,
                t,
                P.Value,
                adj.P.Val)
  

write.table(x=TableOut, file=OutputDE, sep='\t', row.names = F, quote=F)
