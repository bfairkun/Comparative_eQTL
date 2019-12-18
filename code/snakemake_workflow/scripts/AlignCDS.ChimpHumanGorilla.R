library(tidyverse)
library(biomaRt)
library(DECIPHER)


### Get human principal transcripts and lengths

Appris <- read.table("../../data/Appris.principal.isoforms.txt", sep='\t', col.names = c("SYMBOL", "Ensembl.gene", "Ensembl.transcript", "transcript", "apris"), stringsAsFactors = F) %>%
  separate(apris, into=c("apris.type", "apris.score"), remove=F, sep=':')

PrincipalHuamnTranscripts <- Appris %>%
  filter(apris.type=="PRINCIPAL") %>%
  distinct(Ensembl.gene, .keep_all = T)

human_mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# searchAttributes(human_mart, "cds")
print("Got human Mart")
AttributesDesired <- c("ensembl_gene_id",
                       "ensembl_transcript_id",
                       "cds_length"
)

HumanTranscriptLengths <- getBM(attributes = AttributesDesired, mart=human_mart)
print("Got first human getBM query")


HumanTranscriptLengths <- HumanTranscriptLengths %>%
  filter(ensembl_transcript_id %in% PrincipalHuamnTranscripts$Ensembl.transcript)


### Get gorilla orthologs and pick principal transcript

gorilla_mart = useMart("ensembl",dataset="ggorilla_gene_ensembl")
AttributesDesired1 <- c("ensembl_gene_id",
                        "ensembl_transcript_id",
                        "hsapiens_homolog_ensembl_gene",
                        "hsapiens_homolog_orthology_type",
                        "ptroglodytes_homolog_orthology_type"
)
AttributesDesire2 <- c("ensembl_transcript_id", "cds_length")

GorillaQuery1 <- getBM(attributes = AttributesDesired1, mart=gorilla_mart)
GorillaQuery2 <- getBM(attributes = AttributesDesire2, mart=gorilla_mart)

print("Got gorilla getBM queries")


#Combine the queries and filter for orthologs and select transcript closest in length to the human prinicipal transcript
GorillaPrinicpalTranscripts <- GorillaQuery2 %>%
  left_join(GorillaQuery1, by="ensembl_transcript_id") %>%
  filter(
    hsapiens_homolog_orthology_type=="ortholog_one2one" &
      ptroglodytes_homolog_orthology_type=="ortholog_one2one") %>%
  inner_join(HumanTranscriptLengths, by=c("hsapiens_homolog_ensembl_gene"="ensembl_gene_id"), suffix=c(".gorilla", ".human")) %>%
  mutate(DifferenceIn.CDS.Length = abs(cds_length.human-cds_length.gorilla)) %>%
  group_by(hsapiens_homolog_ensembl_gene) %>%
  filter(rank(DifferenceIn.CDS.Length, ties.method="first")==1) %>%
  ungroup()

head(GorillaPrinicpalTranscripts$ensembl_transcript_id.gorilla)

GorillSeqs <- getSequence(
  id=GorillaPrinicpalTranscripts$ensembl_transcript_id.gorilla,
  type='ensembl_transcript_id',
  seqType = 'coding',
  mart=gorilla_mart)
print("Got gorilla seq query")


GorillSeqs.vec <- GorillSeqs$coding
names(GorillSeqs.vec) <- GorillSeqs$ensembl_transcript_id

### Repeat for chimp
chimp_mart = useMart("ensembl",dataset="ptroglodytes_gene_ensembl")
AttributesDesired1 <- c("ensembl_gene_id",
                        "ensembl_transcript_id",
                        "hsapiens_homolog_ensembl_gene",
                        "hsapiens_homolog_orthology_type",
                        "ggorilla_homolog_orthology_type"
)
AttributesDesire2 <- c("ensembl_transcript_id", "cds_length")

ChimpQuery1 <- getBM(attributes = AttributesDesired1, mart=chimp_mart)
ChimpQuery2 <- getBM(attributes = AttributesDesire2, mart=chimp_mart)
print("Got chimp getBM queries")


#Combine the queries and filter for orthologs and select transcript closest in length to the human prinicipal transcript
ChimpPrinicpalTranscripts <- ChimpQuery2 %>%
  left_join(ChimpQuery1, by="ensembl_transcript_id") %>%
  filter(
    hsapiens_homolog_orthology_type=="ortholog_one2one" &
      ggorilla_homolog_orthology_type=="ortholog_one2one") %>%
  inner_join(HumanTranscriptLengths, by=c("hsapiens_homolog_ensembl_gene"="ensembl_gene_id"), suffix=c(".chimp", ".human")) %>%
  mutate(DifferenceIn.CDS.Length = abs(cds_length.human-cds_length.chimp)) %>%
  group_by(hsapiens_homolog_ensembl_gene) %>%
  filter(rank(DifferenceIn.CDS.Length, ties.method="first")==1) %>%
  ungroup()

ChimpSeqs <- getSequence(
  id=ChimpPrinicpalTranscripts$ensembl_transcript_id.chimp,
  type='ensembl_transcript_id',
  seqType = 'coding',
  mart=chimp_mart)
print("Got chimp seq query")


ChimpSeqs.vec <- ChimpSeqs$coding
names(ChimpSeqs.vec) <- ChimpSeqs$ensembl_transcript_id

### And get human seqs too
HumanSeqs <-getSequence(
  id=PrincipalHuamnTranscripts$Ensembl.transcript,
  type='ensembl_transcript_id',
  seqType = 'coding',
  mart=human_mart)
print("Got human seq query")


HumanSeqs.vec <- HumanSeqs$coding
names(HumanSeqs.vec) <- HumanSeqs$ensembl_transcript_id

### Merged CDS, align, and write out aligned seqs to tab file
MergedPrincipalTranscripts <- PrincipalHuamnTranscripts %>%
  dplyr::select(Ensembl.gene, Ensembl.transcript) %>%
  inner_join(
    GorillaPrinicpalTranscripts %>%
      dplyr::select(Ensembl.gene=hsapiens_homolog_ensembl_gene ,ensembl_transcript_id.gorilla),
    by="Ensembl.gene"
  ) %>%
  inner_join(
    ChimpPrinicpalTranscripts %>%
      dplyr::select(Ensembl.gene=hsapiens_homolog_ensembl_gene ,ensembl_transcript_id.chimp),
    by="Ensembl.gene"
  )

print("Aligning CDS and writing out...")
for(i in 1:nrow(MergedPrincipalTranscripts)){
  Gene <- MergedPrincipalTranscripts[i,"Ensembl.gene"]
  HumanID <- MergedPrincipalTranscripts[i,"Ensembl.transcript"]
  ChimpID <- MergedPrincipalTranscripts[i,"ensembl_transcript_id.chimp"]
  GorillaID <- MergedPrincipalTranscripts[i,"ensembl_transcript_id.gorilla"]
  cds <- DNAStringSet(c(
    Chimp=as.character(ChimpSeqs.vec[ChimpID]),
    Human=as.character(HumanSeqs.vec[HumanID]),
    Gorilla=as.character(GorillSeqs.vec[GorillaID])
  ))

  CDS <- AlignTranslation(cds, verbose=FALSE, readingFrame = 1)
  CDS.df <- data.frame(Seq=as.character(CDS), row.names = names(CDS))
  head(CDS.df)
  write.table(CDS.df, file = paste0("~/temp/AlignedCDS.", Gene, ".tab"), col.names=F, quote=F, sep='\t')
}

print("wrote file")


