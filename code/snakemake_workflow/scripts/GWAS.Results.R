setwd("~/CurrentProjects/Comparative_eQTL/code/snakemake_workflow/")
library(data.table)
library(tidyverse)
library(qvalue)
Results <- fread("output/result.assoc.txt.gz", sep='\t')

ManhattanPlot <- Results %>%
  mutate(Color= factor(chr %% 2)) %>%
  filter(chr <= 22) %>%
  ggplot(aes(x=ps, y=-log10(p_wald), color=Color)) +
  geom_point() +
  facet_grid(~chr, scales="free_x", space="free_x") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none") +
  theme(panel.spacing.x=unit(0, "lines"))

ManhattanPlot

ggsave(ManhattanPlot,filename = "../../docs/assets/CellType.GWAS.Manhattan.png", height=2.5, width=8)

#Compare with eQTL results. Start with QQPlot of GWAS P-values grouped by eQTL vs non-eQTL
GTEX.eQTLs <- read.delim("../../output/GTEX_v8_eGenes/Heart_Left_Ventricle.v8.egenes.txt.gz", sep="\t")

EqtlSNPs <- GTEX.eQTLs %>%
  filter(qval<0.01) %>% pull(variant_id)


set.seed(0)
QQ.data <- Results %>%
  dplyr::select(beta, Observed.P=p_wald, rs) %>%
  mutate(IsEQTL = rs %in% EqtlSNPs) %>%
  group_by(IsEQTL) %>%
  sample_n(6159) %>%
  mutate(Expected.P=percent_rank(Observed.P)) %>%
  ungroup()

TestResults <- wilcox.test(Observed.P ~ IsEQTL, data=QQ.data, alternative="less")
lb1 = paste0('P==', format.pval(TestResults$p.value, 2))

QQ.GWAS <-  ggplot(QQ.data, aes(x=-log10(Expected.P), y=-log10(Observed.P), color=IsEQTL)) +
  geom_point() +
  geom_abline() +
  annotate("text",x=-Inf,y=Inf, label=lb1, hjust=-0.1, vjust=1.2, parse=TRUE) +
  theme_bw()
QQ.GWAS
ggsave(QQ.GWAS,filename = "../../docs/assets/CellType.GWAS.QQ.png")

Results$q <- qvalue(Results$p_wald)$qvalues

SigLoci <- Results %>%
  filter(q<0.5)
SigLoci$loci <- SigLoci %>%
  dplyr::select(ps) %>%
  dist(method = "euclidean") %>%
  hclust( method = "complete" ) %>%
  cutree(h=250000)

ClusterLoci <- function(ps){
  tryCatch(
    data.frame(Pos=ps) %>%
      dist(method = "euclidean") %>%
      hclust( method = "complete" ) %>%
      cutree(h=250000) %>%
      return(),
    return(rep(NA, length(ps)))
  )
}

# SigLoci %>% dplyr::select(chr, ps, loci) %>% head(200)

SigLoci %>%
  group_by(chr) %>%
  mutate(Loci=ClusterLoci(ps)) %>% dplyr::select(chr, ps, Loci)

SigLoci %>%
  group_by(loci) %>%
  slice(which.min(q)) %>%
  ungroup() %>%
  dplyr::select(chr, ps, loci, q) %>%
  mutate(stop=ps+1, strand=".") %>%
  dplyr::select(chr, ps, stop, loci, q, strand) %>%
  write_delim("../../output/CellProportionGWAS.loci.bed",delim = '\t',col_names = F)

Results %>%
  sample_n(100) %>%
  dplyr::select(chr, ps, rs, q) %>%
  mutate(stop=ps+1, strand=".") %>%
  dplyr::select(chr, ps, stop, rs, q, strand) %>%
  write_delim("../../output/CellProportionGWAS.RandomControlloci.bed",delim = '\t',col_names = F)
