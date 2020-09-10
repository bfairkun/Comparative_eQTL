setwd("~/CurrentProjects/Comparative_eQTL/code/snakemake_workflow/")
library(data.table)
library(tidyverse)

Results <- fread("output/result.assoc.txt.gz", sep='\t')

ManhattanPlot <- Results %>%
  mutate(Color= factor(chr %% 2)) %>%
  filter(chr <22) %>%
  ggplot(aes(x=ps, y=-log10(p_wald), color=Color)) +
  geom_point() +
  facet_grid(~chr, scales="free_x", space="free_x") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")
ggsave(ManhattanPlot,filename = "../../output/CellType.GWAS.Manhattan.png", height=2.5, width=8)

#Compare with eQTL results. Start with QQPlot of GWAS P-values grouped by eQTL vs non-eQTL
GTEX.eQTLs <- read.delim("../../data/GTEX_v8_eGenes/Heart_Left_Ventricle.v8.egenes.txt.gz", sep="\t")

EqtlSNPs <- GTEX.eQTLs %>%
  filter(qval<0.01) %>% pull(variant_id)

QQ.GWAS <- Results %>%
  dplyr::select(beta, Observed.P=p_wald, rs) %>%
  mutate(IsEQTL = rs %in% EqtlSNPs) %>%
  group_by(IsEQTL) %>%
  sample_n(6159) %>%
  mutate(Expected.P=percent_rank(Observed.P)) %>%
  ungroup() %>%
  ggplot(aes(x=-log10(Expected.P), y=-log10(Observed.P), color=IsEQTL)) +
  geom_point() +
  geom_abline() +
  theme_bw()

ggsave(QQ.GWAS,filename = "../../output/CellType.GWAS.QQ.png")
