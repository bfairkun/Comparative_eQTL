library(tidyverse)
library(data.table)

# setwd("/project2/gilad/bjf79_project1/projects/Comparative_eQTL/code/snakemake_workflow/")
Permuted <-fread("eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.Permuted.txt", drop=c(1:4,6,7))
Real <-fread("eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt", drop=c(1:4,6,7))

# TestPermuted <- data.frame(pvalue=sample(Permuted$pvalue, 5000))
# TestReal <- data.frame(pvalue=sample(Real$pvalue, 5000))

TestFig <- bind_rows(list(Permuted=Permuted, Real=Real), .id = 'source') %>%
  group_by(source) %>%
  mutate(Expected.Pvalue = percent_rank(pvalue)) %>%
  ggplot(aes(y=-log10(pvalue), x=-log10(Expected.Pvalue), color=source)) +
  geom_point() +
  geom_abline() +
  xlab("-log10(Theoretical-Pvalues)") +
  ylab("-log10(Observed-Pvalues)") +
  scale_colour_manual(values=setNames(c("black", "red"), c("Permuted","Real"))) +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  theme(legend.key.width=unit(0.2, "cm"),
        legend.direction = "vertical")
ggsave("../../output/QQPlot.png", TestFig, height=3.7, width=3.7)
