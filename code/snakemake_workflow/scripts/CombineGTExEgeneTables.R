# Script to combine all of GTEx v7 eGene files for all tissues and combine them into a single table.

#setwd("/project/gilad/bjf79/projects/Comparative_eQTL/code/snakemake_workflow/")

files <- list.files(path="scratch/GTEx_Analysis_v7_eQTL", pattern="*.v7.egenes.txt.gz", full.names=TRUE, recursive=FALSE)
Tissues<-gsub("scratch/GTEx_Analysis_v7_eQTL/(.+).v7.egenes.txt.gz", "\\1", files, perl=T)


Merged.eGenes <- data.frame()
for (i in seq_along(files)){
    CurrentTissue <- read.table(files[i], sep='\t', header=T) %>%
      select(gene_id, qval) %>%
      mutate(tissue=Tissues[i])
    Merged.eGenes<-rbind(Merged.eGenes,CurrentTissue)
}

Out <- spread(Merged.eGenes, tissue, qval)
write.table(Out,file = "../../data/AllGTExTissues.egenes.txt",quote=F, sep='\t', row.names = F)

