library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
VepResultFileIn <- args[1]
PnPsFileOut <- args[2]

# VepResultFileIn <- "vep/MergedResults/chimp/ChimpsThisStudy/Merged.txt"
VepResult.df <- fread(VepResultFileIn, select=c(1:7), col.names = c("snp", "pos", "allele", "gene", "transcript", "type", "effect"))


# VepResult.df %>%
#   mutate(MyEffect = case_when(
#     grepl("synonymous",effect) ~ "S",
#     !grepl("synonymous",effect) ~ "N"
#   )) %>%
#   group_by(gene) %>%
#   count(MyEffect) %>%
#   spread(MyEffect, n, fill=0 ) %>%
#   ggplot(aes(x=S, y=N)) +
#   geom_point(alpha=0.1)

VepResult.df %>%
  mutate(MyEffect = case_when(
    grepl("synonymous",effect) ~ "Ps",
    !grepl("synonymous",effect) ~ "Pn"
  )) %>%
  group_by(gene) %>%
  count(MyEffect) %>%
  spread(MyEffect, n, fill=0 ) %>%
  write.table(file=PnPsFileOut, quote=F, sep='\t', row.names = F)