library(TED)
library(tidyverse)
library(biomaRt)
library(Seurat)
source("../code/CustomFunctions.R")



load(system.file("extdata", "tcga.gbm.example.rdata", package="TED"))

ref.norm.4K <- ref.norm %>% head(100)
tcga.tumor.pc.NOchrY.4K <- tcga.tumor.pc.NOchrY %>% head(100)

tcga.ted <- run.Ted (ref.dat = t(ref.norm.4K), X=t(tcga.tumor.pc.NOchrY.4K),
                     tum.key="tumor", input.type="GEP", n.cores=2, outlier.cut=0.01, pdf.name="tcga.tumor")

dim(tcga.ted$res$first.gibbs.res$ theta.merged)
dim(tcga.ted$res$first.gibbs.res$ Znkg)

Heart.seur = readRDS("../big_data/TabMuris_heart_seurat.rds")


#read in table of overdispersion, dispersion, mean expresion estimates
Dispersion <- read.delim('../output/OverdispersionEstimatesFromChimp.txt')

Dispersion$MeanDispersion <- (Dispersion$Chimp.Residual + Dispersion$Human.Residual)/2

CountTableChimpFile <- '../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz'
CountTableHumanFile <- '../output/PowerAnalysisFullCountTable.Human.subread.txt.gz'
OutputDE <- '../output/Final/TableS2.tab'
DropFileName <- '../data/DE_SamplesToDrop.txt'

DropFile <- read.delim(DropFileName, sep='\t', col.names = c("Sample", "Species"), stringsAsFactors = F)
HumanSamplesToDrop <- DropFile %>% filter(Species=="Human") %>% pull(Sample)
ChimpSamplesToDrop <- DropFile %>% filter(Species=="Chimp") %>% pull(Sample)
DE.results <- read.delim(OutputDE, sep='\t', stringsAsFactors = F)
GeneListForOverdispersionCalculation <- DE.results$Ensembl_geneID

CountTables <- GetCountTables(CountTableChimpFile,
                              CountTableHumanFile,
                              0, GeneListForOverdispersionCalculation, ChimpSampleDrop=ChimpSamplesToDrop, HumanSampleDrop = HumanSamplesToDrop)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MouseGenes = rownames(Heart.seur)
genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = MouseGenes ,mart = mouse, attributesL = c("ensembl_gene_id", "chromosome_name", "start_position", "mmusculus_homolog_orthology_type"), martL = human, uniqueRows=T)

# Get list of one to one orthologs
one2one_HumanMouseOrthologs <- genes %>%
  filter(Mouse.homology.type=="ortholog_one2one") %>%
  dplyr::select(MGI.symbol, Gene.stable.ID)

CellTypeSpecificity <- read.table("../output/TissueSpecificity/CellTypeSpecificity.TabulaMurisHeart.tau.log.txt", col.names=c("Gene.stable.ID", "tau"), sep='\t')


CombinedTable <- one2one_HumanMouseOrthologs %>%
  inner_join(CellTypeSpecificity, by="Gene.stable.ID") %>%
  inner_join(Dispersion, by=c("Gene.stable.ID"="gene")) %>%
  drop_na()


#I want to grab the top overdispersed genes, but only those that are reliably expressed in at least 5 cells, and have an ensembl homolog match to human
scCountTable.Filtered <- Heart.seur@assays$RNA[CombinedTable$MGI.symbol,] %>% as.matrix() %>% as.data.frame() %>% head(1000)
CellTypeVector <- Heart.seur@meta.data$Cell.ontology.class
Bulk.CountTable <- cbind(CountTables$Chimp$Counts, CountTables$Human$Counts) %>%
  rownames_to_column() %>%
  right_join(one2one_HumanMouseOrthologs, by=c("rowname"="Gene.stable.ID")) %>%
  filter(MGI.symbol %in% CombinedTable$MGI.symbol) %>%
  dplyr::select(-rowname) %>%
  distinct(MGI.symbol, .keep_all = T) %>%
  column_to_rownames("MGI.symbol") %>% head(1000)



Ted.results <- run.Ted(t(scCountTable.Filtered), t(Bulk.CountTable), pheno.labels = CellTypeVector, input.type = "scRNA", n.cores = 1)
Ted.results$res$first.gibbs.res$Znkg[1,1,] %>% length()
Ted.results$para$X %>% head()
