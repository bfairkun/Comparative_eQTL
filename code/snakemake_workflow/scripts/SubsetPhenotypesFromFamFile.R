args = commandArgs(trailingOnly=TRUE)

# A plink .fam file with N phenotypes, start at column 6 (n=1), to column 6+N (n=N)
FamFileIn <- args[1]
# FamFileIn <- '../eQTL_mapping/plink/ForAssociationTesting.fam'


# File with phenotype numbers (n) to subset.
PhenotypeColumnsFileIn <- args[2]
# PhenotypeColumnsFileIn <- '../scratch/N.txt'

# Output File
OutputFile <- args[3]
# OutputFile <- "testout.txt"

fam.df <- read.table(FamFileIn, sep='\t', stringsAsFactors = F)
phenotypes.to.grab <-  read.table(PhenotypeColumnsFileIn, sep='\t', stringsAsFactors = F, col.names = c('n', 'phenotype_name'))
columns.to.grab <- phenotypes.to.grab$n + 5

MatrixOut <- fam.df[,columns.to.grab]
colnames(MatrixOut) <- phenotypes.to.grab$phenotype_name
rownames(MatrixOut) <- fam.df$V2
write.table(t(MatrixOut), file=OutputFile, quote=F, sep='\t', col.names=NA)
