library(tidyverse)

setwd("/project2/gilad/bjf79_project1/projects/Comparative_eQTL/code/snakemake_workflow/")

LeflerHg38Bed <- read.delim("MiscOutput/Lefler.hg38.bed", col.names =c("chrom", "start", "stop", "name")) %>%
  mutate(ChrPos=paste(chrom, start, sep="."),
         rsID=sub("(rs\\d+).+", "\\1", name))
LeflerHg38Names <- read.delim("MiscOutput/Lefler.hg38.list.txt", col.names =c("chrom", "pos", "name"), sep=" ") %>%
  mutate(ChrPos=paste(chrom, pos, sep="."))
LeflerHg38Frq <- read.delim("MiscOutput/Lefler.hg38.frq", skip=1, col.names=c("chrom", "pos", "N_alleles", "N_chr", "AlleleFreq1", "AlleleFreq2"))%>%
  mutate(ChrPos=paste(chrom, pos, sep="."))
HumanJoined <- inner_join(LeflerHg38Bed, LeflerHg38Names, by="ChrPos", keep=T) %>%
  inner_join(LeflerHg38Frq, by="ChrPos") %>%
  dplyr::select(Human.rsID=rsID, Human.chrom=chrom.x, Human.pos=pos.x, Human.snpName=name.y, AlleleFreq1, AlleleFreq2) %>%
  separate(AlleleFreq1, into=c("Human.RefAllele", "RefAlleleFrq"), sep=":") %>%
  separate(AlleleFreq2, into=c("Human.AltAllele", "Human.AltAlleleFrq"), sep=":") %>%
  dplyr::select(-RefAlleleFrq)



LeflerPanTro5Bed <- read.delim("MiscOutput/Lefler.PanTro5.bed", col.names =c("chrom", "start", "stop", "name")) %>%
  mutate(ChrPos=paste(chrom, start, sep="."),
         rsID=sub("(rs\\d+).+", "\\1", name))
LeflerPanTro5Names <- read.delim("MiscOutput/Lefler.PanTro5.list.txt", col.names =c("chrom", "pos", "name"), sep=" ") %>%
  mutate(ChrPos=paste(chrom, pos, sep="."))
LeflerPanTro5Frq <- read.delim("MiscOutput/Lefler.PanTro5.frq", skip=1, col.names=c("chrom", "pos", "N_alleles", "N_chr", "AlleleFreq1", "AlleleFreq2"))%>%
  mutate(ChrPos=paste(chrom, pos, sep="."))

ChimpJoined <- inner_join(LeflerPanTro5Bed, LeflerPanTro5Names, by="ChrPos", keep=T) %>%
  inner_join(LeflerPanTro5Frq, by="ChrPos") %>%
  dplyr::select(Human.rsID=rsID, Chimp.chrom=chrom.x, Chimp.pos=pos.x, Chimp.snpName=name.y, AlleleFreq1, AlleleFreq2) %>%
  separate(AlleleFreq1, into=c("Chimp.RefAllele", "RefAlleleFrq"), sep=":") %>%
  separate(AlleleFreq2, into=c("Chimp.AltAllele", "Chimp.AltAlleleFrq"), sep=":") %>%
  dplyr::select(-RefAlleleFrq)

full_join(HumanJoined, ChimpJoined, by="Human.rsID") %>%
  write.table("../../output/LeflerTestedSnps.tsv", sep='\t', quote=F, row.names = F)
