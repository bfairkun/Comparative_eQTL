library(TED)
library(tidyverse)
load(system.file("extdata", "tcga.gbm.example.rdata", package="TED"))

ref.norm.4K <- ref.norm %>% head(100)
tcga.tumor.pc.NOchrY.4K <- tcga.tumor.pc.NOchrY %>% head(100)

tcga.ted <- run.Ted (ref.dat = t(ref.norm.4K), X=t(tcga.tumor.pc.NOchrY.4K),
                     tum.key="tumor", input.type="GEP", n.cores=2, outlier.cut=0.01, pdf.name="tcga.tumor")

dim(tcga.ted$res$first.gibbs.res$ theta.merged)
dim(tcga.ted$res$first.gibbs.res$ Znkg)
