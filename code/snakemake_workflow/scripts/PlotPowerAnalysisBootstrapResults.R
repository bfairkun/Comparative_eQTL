library(tidyverse)
library(data.table)
library(pROC)
library(latex2exp)

f <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

Results.f<-"PowerAnalysis/BootstrapRepsMerged.txt.gz"
RealResults.f <- "../../output/Final/TableS2.tab"
ReadSubsetName<-"Full"

args = commandArgs(trailingOnly=TRUE)
Results.f = args[1]
RealResults.f = args[2]
ReadSubsetName = args[3]


Results <- fread(Results.f, sep = '\t', header = F, col.names=c("gene", "log2FC", "P", "P.adjust", "SampleSize", "seed","ReadSubset"))
RealResults <- fread(RealResults.f, sep='\t', header=T)

#Num correct DE genes under two fold change:
Results %>%
  dplyr::select(-c("P")) %>%
  filter(ReadSubset==ReadSubsetName) %>%
  # filter(seed %in% 1:10) %>%
  left_join(
    (RealResults %>%
       dplyr::select(Ensembl_geneID, logFC)),
    suffix=c(".Real", ".Resampled"),
    by=c("gene"="Ensembl_geneID")
  ) %>%
  mutate(FDR.01=P.adjust<0.01,
         FDR.05=P.adjust<0.05,
         FDR.10=P.adjust<0.1) %>%
  gather(key="FDR", value = "TestResults", -c("gene", "P.adjust", "SampleSize", "seed", "ReadSubset", "logFC", "log2FC") ) %>%
  mutate(FDR=recode(FDR, FDR.01="FDR<0.01", FDR.05="FDR<0.05", FDR.10="FDR<0.10")) %>%
  mutate(NewTestResults=(TestResults & abs(logFC)<=1)) %>%
  group_by(FDR, SampleSize, seed) %>%
  summarize(NumDiscoveries=sum(NewTestResults, na.rm=T)) %>%
  group_by(FDR, SampleSize) %>%
  summarize(med=median(NumDiscoveries, na.rm=T),
            se=sd(NumDiscoveries, na.rm=T),
            LQ=quantile(NumDiscoveries, 0.25, na.rm=T),
            UQ=quantile(NumDiscoveries, 0.75, na.rm=T)) %>%
  filter(SampleSize==4)




#Num DE genes
ToPlot <- Results %>%
  dplyr::select(-c("log2FC", "P")) %>%
  filter(ReadSubset==ReadSubsetName) %>%
  mutate(FDR.01=P.adjust<0.01,
         FDR.05=P.adjust<0.05,
         FDR.10=P.adjust<0.1) %>%
  gather(key="FDR", value = "TestResults", -c("gene", "P.adjust", "SampleSize", "seed", "ReadSubset") ) %>%
  mutate(FDR=recode(FDR, FDR.01="FDR<0.01", FDR.05="FDR<0.05", FDR.10="FDR<0.10")) %>%
  group_by(FDR, SampleSize, seed) %>%
  summarize(NumDiscoveries=sum(TestResults))

ToPlotMed <- ToPlot %>%
  group_by(SampleSize, FDR) %>%
  summarize(MedianNumDiscoveries=median(NumDiscoveries),
            IQR=IQR(NumDiscoveries),
            avg=mean(NumDiscoveries, na.rm=T),
            se=sd(NumDiscoveries, na.rm=T),
            LQ=quantile(NumDiscoveries, 0.25, na.rm=T),
            UQ=quantile(NumDiscoveries, 0.75, na.rm=T))


TrueNumDiscoveries <- RealResults %>%
  dplyr::select(P.adjust=adj.P.Val) %>%
  mutate(FDR.01=P.adjust<0.01,
         FDR.05=P.adjust<0.05,
         FDR.10=P.adjust<0.1) %>%
  gather(key="FDR", value = "TestResults", -P.adjust) %>%
  mutate(FDR=recode(FDR, FDR.01="FDR<0.01", FDR.05="FDR<0.05", FDR.10="FDR<0.10")) %>%
  group_by(FDR) %>%
  summarize(NumDiscoveries=sum(TestResults))


ggplot(ToPlot, aes(x=SampleSize, y=NumDiscoveries, group=interaction(SampleSize, FDR), color=FDR)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_line(data=ToPlotMed, aes(x=SampleSize, y=MedianNumDiscoveries, group=FDR)) +
  stat_summary(fun.data = f, geom="boxplot", position="dodge2") +
  geom_hline(data=TrueNumDiscoveries, aes(yintercept=NumDiscoveries, color=FDR), linetype='dashed') +
  scale_x_continuous(limits=c(1,40), breaks=c(2,4,8,12,16,24,32,39)) +
  ylab("Number DE genes") +
  theme_bw() +
  theme(legend.key = element_rect(colour = "transparent", fill = NA),
        legend.justification = c(1,0), legend.position = c(0.9,0),
        legend.background=element_blank(),
        legend.title=element_blank())

ggsave(paste0("../../output/Final/FigS_Power_NumDE_", ReadSubsetName, ".pdf"), width=3, height=3)

# percent false discoveries
FDR.ToPlot <- Results %>%
  dplyr::select(-P) %>%
  filter(ReadSubset==ReadSubsetName) %>%
  left_join(
    (RealResults %>%
       dplyr::select(Ensembl_geneID, Real.P.adjust=adj.P.Val, Real.FC=logFC)),
    suffix=c(".Real", ".Resampled"),
    by=c("gene"="Ensembl_geneID")
  ) %>%
  mutate(FDR.01=P.adjust<0.01,
         FDR.05=P.adjust<0.05,
         FDR.10=P.adjust<0.1,
         Real.FDR.01=Real.P.adjust<0.01) %>%
  dplyr::select(-c("Real.P.adjust", "P.adjust", "ReadSubset")) %>%
  gather(key="FDR", value = "TestResults", -c("gene", "SampleSize", "seed", "Real.FDR.01", "Real.FC", "log2FC")) %>%
  #False discoveries include ones where the effect direction are incorrect
  mutate(FalseDiscovery=(TestResults & !Real.FDR.01) | (TestResults & (sign(Real.FC) != sign(log2FC) ))) %>%
  # mutate(FalseDiscovery=(TestResults & !Real.FDR.01)) %>%a
  drop_na() %>%
  group_by(FDR, SampleSize, seed) %>%
  summarize(FalseDiscoveryRate=sum(FalseDiscovery)/sum(TestResults))

FDR.ToPlotMed <- FDR.ToPlot %>%
  group_by(SampleSize, FDR) %>%
  summarize(MedianFDR=median(FalseDiscoveryRate),
            LQ=quantile(FalseDiscoveryRate, 0.25, na.rm=T),
            UQ=quantile(FalseDiscoveryRate, 0.75, na.rm=T),
            avg=mean(FalseDiscoveryRate, na.rm=T),
            se=sd(FalseDiscoveryRate, na.rm=T))


ggplot(FDR.ToPlot, aes(x=SampleSize, y=FalseDiscoveryRate, group=interaction(SampleSize, FDR), color=FDR)) +
  # geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.data = f, geom="boxplot", position="dodge2") +
  scale_x_continuous(limits=c(1,40), breaks=c(2,4,8,12,16,24,32,39)) +
  scale_y_continuous(limits=c(0,0.3)) +
  # geom_line(data=FDR.ToPlotMed, aes(x=SampleSize, y=MedianFDR, group=FDR)) +
  ylab("Empricial FDR estimate") +
  theme_bw() +
  theme(legend.key = element_rect(colour = "transparent", fill = NA),
        legend.justification = c(1,1), legend.position = c(0.9,0.9),
        legend.background=element_blank(),
        legend.title=element_blank())
ggsave(paste0("../../output/Final/FigS_Power_FDREstimate_", ReadSubsetName, ".pdf"), width=3, height=3)



#Plot of distribution of true effect sizes of true positives at different FDR.
## First get median adj.P.value across bootstrap reps
Results %>%
  filter(ReadSubset==ReadSubsetName) %>%
  group_by(SampleSize, gene) %>%
  drop_na() %>%
  summarize(Median.P.adjust=median(P.adjust)) %>%
  mutate(FDR.01=Median.P.adjust<0.01,
         FDR.05=Median.P.adjust<0.05,
         FDR.10=Median.P.adjust<0.1) %>%
  dplyr::select(-c("Median.P.adjust")) %>%
gather(key="FDR", value = "TestResults", -c("gene", "SampleSize")) %>%
  left_join(RealResults, by=c("gene"="Ensembl_geneID")) %>%
  filter(TestResults==TRUE & adj.P.Val < 0.01) %>%
  ggplot(aes(x=SampleSize, y=abs(logFC), group=interaction(SampleSize, FDR), color=FDR)) +
  stat_summary(fun.data = f, geom="boxplot", position="dodge2") +
  scale_x_continuous(limits=c(1,40), breaks=c(2,4,8,12,16,24,32,39)) +
  ylab(expression(atop("Effect-size of true positive DE genes",'|log'[2]*'(FoldChange)|'))) +
  xlab("Sample Size") +
  theme_bw() +
  theme(legend.key = element_rect(colour = "transparent", fill = NA),
        legend.justification = c(1,1), legend.position = c(0.9,0.9),
        legend.background=element_blank(),
        legend.title=element_blank())
ggsave(paste0("../../output/Final/FigS_Power_EffectSize_", ReadSubsetName, ".pdf"), width=3, height=3)


ReadSubsetName<-"Full"

#Plot median difference in effect size estimates of DE genes estimated using full set versus subset.
#I expect a larger bias (Winner's curse) in smaller sample sizes
WinnersCursePlotData <- Results %>%
  filter(ReadSubset==ReadSubsetName) %>%
  drop_na() %>%
  dplyr::select("gene", "log2FC", "P.adjust", "SampleSize", "seed") %>%
  mutate(FDR.01=P.adjust<0.01,
         FDR.05=P.adjust<0.05,
         FDR.10=P.adjust<0.1) %>%
  gather(key="FDR", value = "TestResults", -c("gene", "SampleSize", "log2FC", "seed", "P.adjust")) %>%
  filter(TestResults==T) %>%
  left_join((RealResults %>%
               dplyr::select(Ensembl_geneID, TrueLog2FC=logFC, True.adj.P.val=adj.P.Val)),
            by=c("gene"="Ensembl_geneID")) %>%
  mutate(DifferenceInEffectSizeEstimate=log2FC*sign(TrueLog2FC)-TrueLog2FC*sign(TrueLog2FC)) %>%
  group_by(SampleSize, FDR, seed) %>%
  summarize(MedianDifference = median(DifferenceInEffectSizeEstimate, na.rm = T))
ggplot(WinnersCursePlotData, aes(x=SampleSize, y=MedianDifference, group=interaction(SampleSize, FDR), color=FDR)) +
  stat_summary(fun.data = f, geom="boxplot", position="dodge2") +
  ylim(c(0,1)) +
  ylab("Median magnitude of effect size over-estimates") +
  theme_bw() +
  theme(legend.key = element_rect(colour = "transparent", fill = NA),
        legend.justification = c(1,1), legend.position = c(0.9,0.9),
        legend.background=element_blank(),
        legend.title=element_blank())
ggsave(paste0("../../output/Final/FigS_Power_WinnersCurse_", ReadSubsetName, ".pdf"), width=3, height=3)





#ROC curves
RocInfo <- Results %>%
  dplyr::select(-c("log2FC", "P.adjust")) %>%
  filter(ReadSubset==ReadSubsetName) %>%
  # filter(seed %in% 1:10) %>%
  left_join(
    (RealResults %>%
       dplyr::select(Ensembl_geneID, Real.P=P.Value, Real.P.adjust=adj.P.Val)),
    suffix=c(".Real", ".Resampled"),
    by=c("gene"="Ensembl_geneID")
  ) %>%
  mutate(TrueResponse=Real.P.adjust<0.01) %>%
  dplyr::select(-Real.P.adjust) %>%
  unite(SampleSize_seed, SampleSize, seed) %>%
  drop_na()



#Reformat data to make it easier to iterate through regions and calculate smoothed densities
group_names <- RocInfo %>%
  group_keys(SampleSize_seed) %>% pull(1)
ListOfDf <- RocInfo %>%
  group_split(SampleSize_seed) %>%
  set_names(group_names)

#Calculate roc objects.
Roc.results <- lapply(ListOfDf,function(i) {roc(response=i$TrueResponse, predictor=i$P, direction=">", quiet = T)})

#Get coordinates of specificity, sensitivity from roc objects
Roc.coords <- lapply(Roc.results, function(i) {coords(i, x=c(1:10 %o% 10^(-50:0)), transpose=F, ret=c("sensitivity", "specificity", "threshold", "fdr"))})

#Convert list to df and get median, and 5, 95 percentile
Roc.coords.df <- bind_rows(Roc.coords, .id = 'SampleSize_seed') %>%
  separate(SampleSize_seed, into=c("SampleSize", "seed"), convert=T, sep="_") %>%
  group_by(SampleSize, threshold) %>%
  summarize(sensitivity.lower=quantile(sensitivity, probs=0.05),
            sensitivity.med=median(sensitivity),
            sensitivity.upper=quantile(sensitivity, probs=0.95),
            specificity.med=median(specificity),
            fdr.med=median(fdr, na.rm = T))

Roc.coords.df %>% filter(!SampleSize %in% c(12, 32)) %>%
ggplot(aes(x=1-specificity.med, y=sensitivity.med, color=as.factor(SampleSize))) +
  geom_abline() +
  geom_line() +
  xlab("false positive rate (1-specificity)") +
  ylab("true positive rate (sensitivity)") +
  scale_colour_discrete(name  ="Sample size") +
  labs(color = "Sample size") +
  geom_ribbon(aes(ymin=sensitivity.lower, ymax=sensitivity.upper), linetype="dashed", alpha=0.0) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.key = element_rect(colour = "transparent", fill = NA),
        legend.justification = c(1,0), legend.position = c(0.95,0.05),
        legend.background=element_blank()) +
  guides(color=guide_legend(ncol=2))
ggsave(paste0("../../output/Final/FigS_Power_ROC_", ReadSubsetName, ".pdf"), width=3, height=3)


bind_rows(Roc.coords, .id = 'SampleSize_seed') %>%
  separate(SampleSize_seed, into=c("SampleSize", "seed"), convert=T, sep="_") %>%
  filter(!SampleSize %in% c(12, 32)) %>%
  filter(seed %in% c(1:10)) %>%
  ggplot(aes(x=1-specificity, y=sensitivity, group=interaction(as.factor(SampleSize), seed), color=as.factor(SampleSize))) +
  geom_abline() +
  geom_line() +
  scale_color_brewer(palette = "Set2") +
  labs(color = "Sample size") +
  xlab("false positive rate (1-specificity)") +
  ylab("true positive rate (sensitivity)") +
  theme_bw() +
  theme(legend.key = element_rect(colour = "transparent", fill = NA),
        legend.justification = c(1,0), legend.position = c(0.95,0.05),
        legend.background=element_blank()) +
  guides(color=guide_legend(ncol=2))
ggsave(paste0("../../output/Final/FigS_Power_10ROCs_", ReadSubsetName, ".pdf"), width=3, height=3)


#ROC like plot with fdr on x-axis
Roc.coords.df %>% filter(!SampleSize %in% c(12, 32)) %>%
  filter(fdr.med <0.5) %>%
  ggplot(aes(x=fdr.med, y=sensitivity.med, color=as.factor(SampleSize))) +
  geom_line() +
  xlab("FDR") +
  ylab("true positive rate (sensitivity)") +
  scale_colour_discrete(name  ="Sample size") +
  labs(color = "Sample size") +
  geom_ribbon(aes(ymin=sensitivity.lower, ymax=sensitivity.upper), linetype="dashed", alpha=0.0) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.key = element_rect(colour = "transparent", fill = NA),
        legend.justification = c(1,0), legend.position = c(0.99,0.01),
        legend.background=element_blank()) +
  guides(color=guide_legend(ncol=2))
ggsave(paste0("../../output/Final/FigS_Power_SensitivityFDR_", ReadSubsetName, ".pdf"), width=3, height=3)




