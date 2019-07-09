library(tidyverse)

LDDecay <- read.table(gzfile("~/CurrentProjects/Comparative_eQTL/code/snakemake_workflow/scratch/LDDecay.txt.stat.gz"), header=T, comment.char="@")
ggplot(LDDecay, aes(y=Mean_r.2, x=X.Dist)) +
  geom_line()
