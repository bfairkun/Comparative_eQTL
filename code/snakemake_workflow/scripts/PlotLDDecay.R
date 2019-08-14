library(tidyverse)

LDDecay <- read.table(gzfile("~/CurrentProjects/Comparative_eQTL/code/snakemake_workflow/scratch/LDDecay.txt.stat.gz"), header=T, comment.char="@")
ggplot(LDDecay, aes(y=Mean_r.2, x=X.Dist)) +
  geom_line()

LDDecayYRI <- read.table(gzfile("~/CurrentProjects/Comparative_eQTL/code/snakemake_workflow/scratch/YRI.LDDecay.stat.gz"), header=T, comment.char="@")
ggplot(LDDecayYRI, aes(y=Mean_r.2, x=X.Dist)) +
  geom_line()


# LDDecayChr1 <- read.table(gzfile("~/CurrentProjects/Comparative_eQTL/code/snakemake_workflow/scratch/LDDecay.ChimpChr1.txt.stat.gz"), header=T, comment.char="@")
ggplot(LDDecay, aes(y=Mean_r.2, x=X.Dist)) +
  geom_line() +
  geom_line(data=LDDecayYRI, color="red") +
  xlim(c(1,100000)) +
  xlab("SNP Distnace") +
  ylab("Mean R-squared") +
  theme_bw()

