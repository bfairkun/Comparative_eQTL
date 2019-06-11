library("clusterProfiler")
library("org.Hs.eg.db")

bitr(c("ENSG00000196611", "ENSG00000093009"),
     fromType = "ENSEMBL",
     toType = c("ENTREZID", "SYMBOL"),
     OrgDb = org.Hs.eg.db)
