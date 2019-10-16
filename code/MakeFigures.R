library(tidyverse)
library(cowplot)

source("code/CustomFunctions.R")

A<-TsvToCombinedEgenes(Chimp.tsv = "output/ChimpEgenes.eigenMT.txt.gz", Human.tsv = "output/GTEX_renalysis/SampleSize_80.txt.gz", SysToID.tsv = "data/Biomart_export.Hsap.Ptro.orthologs.txt.gz")
B<-TsvToCombinedEgenes(Chimp.tsv = "output/ChimpEgenes.eigenMT.txt.gz", Human.tsv = "output/GTEX_renalysis/SampleSize_60.txt.gz", SysToID.tsv = "data/Biomart_export.Hsap.Ptro.orthologs.txt.gz")

Test <- Plot.PercentNonIdentity.byGroup(A)
IdentityPlot <- ggdraw(Test$plot) + draw_grob(tableGrob(Test$PvalTable, rows=NULL, theme=ttheme_default(colhead=list(fg_params = list(parse=TRUE)))), x=0.7, y=0.15, width=0.2, height=0.3)


Test.dNdS <- Plot.dNdS.byGroup(A)
dNdSPlot <- ggdraw(Test.dNdS$plot) + draw_grob(tableGrob(Test.dNdS$PvalTable, rows=NULL, theme=ttheme_default(colhead=list(fg_params = list(parse=TRUE)))), x=0.7, y=0.15, width=0.2, height=0.3)
ToSave <- plot_grid(IdentityPlot, dNdSPlot, labels = c('A', 'B'), scale=c(1,1))
# ggsave(filename="~/Desktop/Test.pdf", ToSave, height=4, width=8)

OddsRatiosByHumanSampleSize.df <- data.frame()
for (filename in Sys.glob("./output/GTEX_renalysis/SampleSize_*.txt.gz")){
  E.genes.df <- TsvToCombinedEgenes(Chimp.tsv = "output/ChimpEgenes.eigenMT.txt.gz", Human.tsv = filename, SysToID.tsv = "data/Biomart_export.Hsap.Ptro.orthologs.txt.gz")
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
    SampleSize=gsub("./output/GTEX_renalysis/SampleSize_(\\d+).txt.gz", "\\1",filename, perl=T) %>% as.numeric()
  )
  OddsRatiosByHumanSampleSize.df <- rbind(OddsRatiosByHumanSampleSize.df, E.genes.OddsRatios)
}
ggplot(OddsRatiosByHumanSampleSize.df, aes(x=NumHumanEgenes, y=OR))+
  geom_point(aes(size=SampleSize)) +
  geom_hline(yintercept=1, linetype="dashed") +
  ylab("Odds Ratio") +
  xlab("Number eGenes discovered in human") +
  theme_bw()
