library(tidyverse)
library(cowplot)
library('latex2exp')
setwd("/Users/benfair/Documents/GiladLiProjects/Repos/Comparative_eQTL/code/snakemake_workflow/")

source("../CustomFunctions.R")

# Fig1 --------------------------------------------------------------------

#A: example plot of highly and lowly dispersed genes
#get example genes
ParameterEstimates <- read.table("../../output/OverdispersionEstimatesFromChimp_NoVirusChallangedIndividuals.txt.gz", sep='\t', header=T, stringsAsFactors = F)
HighAndLowDispersionGenes <- ParameterEstimates %>%
  filter(Chimp.Mean.Log2RPKM >= 2 & Chimp.Mean.Log2RPKM <= 3) %>%
  arrange(Chimp.Overdispersion) %>%
  filter(row_number()==1 | row_number()==n()) %>% pull(gene)

#Read count table
HumanSamplesToDrop <- c(c("SRR1507229","SRR603918", "SRR1478149", "SRR598509", "SRR613186"), c("SRR1489693", "SRR598148", "59167", "SRR1478900", "SRR1474730", "61317"))
ChimpSamplesToDrop <- c("Little_R")
OtherMetadata <- as.data.frame(read_excel("../../data/Metadata.xlsx"))
VirusChallengedChimps <- OtherMetadata %>% filter(grepl("V+",Viral.status)) %>% pull(Individual.ID)
CountTableChimpFile <- '../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz'
CountTableHumanFile <- '../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz'

EgenesTested <- TsvToCombinedEgenes(Chimp.tsv = "../../output/ChimpEgenes.eigenMT.txt.gz", Human.tsv = "../../data/GTEX_v8_eGenes/Heart_Left_Ventricle.v8.egenes.txt.gz", SysToID.tsv = "../../data/Biomart_export.Hsap.Ptro.orthologs.txt.gz", HumanTsvType = "GTEx")
GeneListTestedForQTLs <- EgenesTested$H.gene
GenesInDESet <- read.table(gzfile(CountTableChimpFile), header=T, check.names=FALSE, skip=1)$Geneid
GeneListForOverdispersionCalculation <- intersect(as.character(GenesInDESet),GeneListTestedForQTLs)

CountTables <- GetCountTables(CountTableChimpFile,
                              CountTableHumanFile,
                              0, GeneListForOverdispersionCalculation, ChimpSampleDrop=ChimpSamplesToDrop, HumanSampleDrop = HumanSamplesToDrop)

Fig1ATable <- ParameterEstimates %>% filter(gene %in% HighAndLowDispersionGenes) %>%
  dplyr::select(gene, Chimp.Mean.Expression, Chimp.Overdispersion)
colnames(Fig1ATable) <- c('gene', '"mean," ~ hat(mu)', '"overdispersion," ~ hat(phi)')
Fig1A <- CountTables$Chimp$log2RPKM[HighAndLowDispersionGenes,] %>% t() %>% as.data.frame() %>% gather() %>%
  dplyr::select(gene=key, value) %>%
  ggplot(aes(x=gene, y=value)) +
  geom_boxplot(outlier.shape=NA) +
  ylab("Expression\nlog(RPKM)") +
  geom_jitter(position=position_jitter(width=.1, height=0))
Fig1A_withTable <- ggdraw(Fig1A) +
  draw_grob(tableGrob(Fig1ATable, rows=NULL, theme=ttheme_default(colhead=list(fg_params = list(parse=TRUE)))), x=0.4, y=0.7, width=0.2, height=0.3)

CountTables$Chimp$log2RPKM[HighAndLowDispersionGenes,] %>% t() %>% as.data.frame() %>% gather() %>%
  dplyr::select(gene=key, value) %>%
  ggplot(aes(x=gene, y=value)) +
  geom_boxplot(outlier.shape=NA) +
  ylab("Expression\nlog2(RPKM)") +
  ylim(c(0,10)) +
  geom_jitter(position=position_jitter(width=.1, height=0))


#overdispersion correlation
R<-cor(log(ParameterEstimates$Chimp.Overdispersion),log(ParameterEstimates$Human.Overdispersion), use="complete.obs")
lb1 <- paste("~R^2==~", round(R**2,2))
Fig1B <- ggplot(ParameterEstimates, aes(x=Chimp.Overdispersion, y=Human.Overdispersion)) +
  geom_point(alpha=0.05) +
  scale_x_continuous(trans='log10', limits=c(1E-2,10), name=expression(paste("Chimp overdispersion, ", "1","/", hat(phi)))) +
  scale_y_continuous(trans='log10', limits=c(1E-2,10), name=expression(paste("Human overdispersion, ", "1","/", hat(phi)))) +
  annotate("text",x=10,y=.01, label=lb1, hjust=1, vjust=-1, parse=TRUE) +
  theme_bw()

#Even when regressing out the mean effect
LoessPlot <- ggplot(ParameterEstimates, aes(x=Chimp.Mean.Expression, y=Chimp.Overdispersion)) +
  geom_point(alpha=0.1, size=0.5) +
  scale_x_continuous(name=TeX('$\\hat{\\mu}$'), limits=c(-25,-12.5)) +
  scale_y_continuous(trans="log10", name=expression(paste("log(", hat(phi), ")")), limits=c(0.01,10)) +
  geom_smooth(method=loess, show.legend = FALSE, se=T, method.args=list(degree=1, family="gaussian")) +
  # geom_smooth(method=lm, show.legend = FALSE, se=T) +
  theme_bw()
# LoessPlot
ResidualDemoPlot <- LoessPlot +
  geom_segment(aes(x=-21.125,xend=-21.125, y=0.35, yend=2.5),
               lineend = "round", linejoin = "round", color="brown",
               size = 0.7, arrow = arrow(length = unit(0.1, "inches"))
  ) +
  # geom_density_2d() +
  annotate("text",x=-21.125,y=2.5, label=TeX("$\\epsilon_j$"), hjust=0, vjust=0, color="brown", size=10) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = 'transparent', colour = 'transparent'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ResidualDemoPlot

R<-cor(log(ParameterEstimates$Chimp.Overdispersion),log(ParameterEstimates$Human.Overdispersion), use="complete.obs")
lb1 <- paste("~R^2==~", round(R**2,2))
CorrelationAfterRegressingOut <- ggplot(ParameterEstimates, aes(x=Chimp.Residual, y=Human.Residual)) +
  geom_point(alpha=0.05) +
  scale_x_continuous(limits=c(0.05,50), trans="log10", name=TeX("Chimp residual, $\\epsilon_j$")) +
  scale_y_continuous(limits=c(0.05,50), trans="log10", name=TeX("Human residual, $\\epsilon_j$")) +
  annotate("text",x=50,y=0.05, label=lb1, hjust=1, vjust=-1, parse=TRUE) +
  theme_bw()

ResidualsCorrelateWithInset <- ggdraw(CorrelationAfterRegressingOut) +
  draw_plot(ResidualDemoPlot, .1, .5, .4, .5)


#heatmap with gtex




# FigX --------------------------------------------------------------------


A<-TsvToCombinedEgenes(Chimp.tsv = "../../output/ChimpEgenes.eigenMT.txt.gz", Human.tsv = "../../output/GTEX_renalysis/SampleSize_80.txt.gz", SysToID.tsv = "../../data/Biomart_export.Hsap.Ptro.orthologs.txt.gz")
B<-TsvToCombinedEgenes(Chimp.tsv = "../../output/ChimpEgenes.eigenMT.txt.gz", Human.tsv = "../../output/GTEX_renalysis/SampleSize_60.txt.gz", SysToID.tsv = "../../data/Biomart_export.Hsap.Ptro.orthologs.txt.gz")

Test <- Plot.PercentNonIdentity.byGroup(A) + ylab("Percent non-identity between\nchimp and human")
IdentityPlot <- ggdraw(Test$plot) + draw_grob(tableGrob(Test$PvalTable, rows=NULL, theme=ttheme_default(colhead=list(fg_params = list(parse=TRUE)))), x=0.7, y=0.3, width=0.2, height=0.3)


Test.dNdS <- Plot.dNdS.byGroup(A)
dNdSPlot <- ggdraw(Test.dNdS$plot) + draw_grob(tableGrob(Test.dNdS$PvalTable, rows=NULL, theme=ttheme_default(colhead=list(fg_params = list(parse=TRUE)))), x=0.7, y=0.3, width=0.2, height=0.3)
ToSave <- plot_grid(IdentityPlot, dNdSPlot, labels = c('A', 'B'), scale=c(1,1))
ggsave(filename="~/Desktop/Test.pdf", ToSave, height=6, width=12)

OddsRatiosByHumanSampleSize.df <- data.frame()
for (filename in Sys.glob("../../output/GTEX_renalysis/SampleSize_*.txt.gz")){
  E.genes.df <- TsvToCombinedEgenes(Chimp.tsv = "../../output/ChimpEgenes.eigenMT.txt.gz", Human.tsv = filename, SysToID.tsv = "../../data/Biomart_export.Hsap.Ptro.orthologs.txt.gz")
  ContigencyTable <- matrix( c( sum(E.genes.df$H.FDR<=0.1 & E.genes.df$C.FDR<=0.1),
                                sum(E.genes.df$H.FDR<=0.1 & E.genes.df$C.FDR>0.1),
                                sum(E.genes.df$H.FDR>0.1 & E.genes.df$C.FDR<=0.1),
                                sum(E.genes.df$H.FDR>0.1 & E.genes.df$C.FDR>0.1)),
                             nrow = 2)

  rownames(ContigencyTable) <- c("Chimp eGene", "Not Chimp eGene")
  colnames(ContigencyTable) <- c("Human eGene", "Not human eGene")
  FisherResults <- (fisher.test(ContigencyTable)
)
  E.genes.OddsRatios <- data.frame(
    OR=FisherResults$estimate,
    Upper=FisherResults$conf.int[[2]],
    Lower=FisherResults$conf.int[[1]],
    P=fisher.test(ContigencyTable, alternative="greater")$p.value,
    NumHumanEgenes=sum(E.genes.df$H.FDR<=0.1),
    SampleSize=gsub("../../output/GTEX_renalysis/SampleSize_(\\d+).txt.gz", "\\1",filename, perl=T) %>% as.numeric()
  )
  OddsRatiosByHumanSampleSize.df <- rbind(OddsRatiosByHumanSampleSize.df, E.genes.OddsRatios)
}
ggplot(OddsRatiosByHumanSampleSize.df, aes(x=NumHumanEgenes, y=OR))+
  geom_point(aes(size=SampleSize)) +
  geom_hline(yintercept=1, linetype="dashed") +
  ylab("Odds Ratio") +
  xlab("Number eGenes discovered in human") +
  theme_bw()
