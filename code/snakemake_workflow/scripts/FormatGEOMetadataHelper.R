library(tidyverse)
library(readxl)

md5sums <- read.delim("DataForSubmission/GEO/FastqChecksums.txt", col.names=c("checksum", "blank", "fastq"), sep = " ", stringsAsFactors = F) %>%
  mutate(basename=basename(fastq)) %>%
  select(basename, checksum)
Samples <- read.delim("../../data/NovelRNASeqFiles.tsv", stringsAsFactors = F) %>%
  mutate(basename=basename(fastq)) %>%
  select(sample, basename) %>%
  group_by(sample) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  spread(id, basename)

write_delim(Samples, "DataForSubmission/GEO/MetadataSamples.tsv", delim = '\t')
write_delim(md5sums, "DataForSubmission/GEO/MetadataChecksums.tsv", delim='\t')

Metadata <- read_excel("../../data/Metadata.xlsx")
Units <- read.delim("units.tsv", stringsAsFactors = F) %>%
  mutate(
    fq1=basename(fq1),
    fq2=basename(fq2)
  )
Metadata %>%
  left_join(Units, by=c("Individual.ID"="sample")) %>%
  dplyr::select(Individual.ID, fq1, fq2) %>%
  write_delim("DataForSubmission/GEO/WGS.SRA.Fqs.tsv", delim = '\t')
