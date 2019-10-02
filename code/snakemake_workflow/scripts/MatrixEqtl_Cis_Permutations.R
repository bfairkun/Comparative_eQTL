# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
SNP_file_name <- args[1]
snps_location_file_name <- args[2]
expression_file_name <- args[3]
gene_location_file_name <- args[4]
covariates_file_name <- args[5]
errorCovariance_file <- args[6]
permutation_matrix_output_filename <- args[7]
Npermutations <- as.numeric(args[8])
InitialSeed <- as.numeric(args[9])
cisDistance <- args[10]

# setwd("/project2/gilad/bjf79_project1/projects/Comparative_eQTL/")
# SNP_file_name <- "code/snakemake_workflow/scratch/Test.snps"
# snps_location_file_name <- "code/snakemake_workflow/scratch/Test.snploc"
# expression_file_name <- "code/snakemake_workflow/eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt"
# gene_location_file_name <- "code/snakemake_workflow/eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt"
# covariates_file_name <- "output/Covariates/0GenotypePCs_and_11RNASeqPCs.covariates"
# errorCovariance_file <- "code/snakemake_workflow/eQTL_mapping/Kinship/GRM.cXX.txt"
# permutation_matrix_output_filename <- "/project2/gilad/bjf79/temp/PvalueMatrix.txt"
# cisDistance<-100000
# InitialSeed <- 1
# Npermutations <- 5


# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS


output_file_name_tra = tempfile();

# Only associations significant at this level will be saved

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance <- as.matrix(read.table(errorCovariance_file,sep='\t'))


# Distance for local gene-SNP pairs
cisDist = as.numeric(cisDistance);
print(paste("Cis distance setting:",cisDist))

## Load genotype data
print("Loading genotype data...")
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data
print("Loading expression data...")
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
print("Loading covariate data...")
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Other initialization stuff
print("Loading other stuff...")
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

### Permute the sample labels for expression ###
# Read in actual data expression matrix and covariates
ActualData.ExpressionMatrix <- read.table(expression_file_name, header=T, row.names = 1, check.names = FALSE )
ActualData.Cov <- read.table(covariates_file_name, header=T, row.names = 1, check.names = FALSE )

#Matrix. Rows are genes. Columns are permutations
PermutatePvalueMatrix <- matrix(data=NA, nrow = nrow(ActualData.ExpressionMatrix), ncol = Npermutations)

### Run permutations
print("Running permutations...")
for (i in 1:Npermutations){
  # Permute column labels for both (using same seed for randomization)
  # technically I am permuting the column data, and preserving the column labels. That way, I do not have to do the same for the (huge) genotype file, which would take more computational time
  set.seed(i+InitialSeed)
  Temp.df <- ActualData.ExpressionMatrix %>% select(sample(colnames(ActualData.ExpressionMatrix), length(colnames(ActualData.ExpressionMatrix))))
  set.seed(i+InitialSeed)
  Temp.df.cov <- ActualData.Cov %>% select(sample(colnames(ActualData.ExpressionMatrix), length(colnames(ActualData.ExpressionMatrix))))
  colnames(Temp.df) <- colnames(ActualData.ExpressionMatrix)
  colnames(Temp.df.cov) <- colnames(ActualData.ExpressionMatrix)

  # Write out permutated expression matrix and reload it
  TempFilepath.ExpressionMatrix <- tempfile("ExpressionMatrix.")
  write.table(Temp.df, file=TempFilepath.ExpressionMatrix, sep='\t', quote=F, col.names =NA)
  gene$LoadFile(TempFilepath.ExpressionMatrix);

  TempFilepath.Covariates <- tempfile("Covariates.")
  write.table(Temp.df.cov, file=TempFilepath.Covariates, sep='\t', quote=F, col.names =NA)
  cvrt$LoadFile(TempFilepath.Covariates);

  #Calculate Pvalues from permutated data
  permuted = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = NULL,
    pvOutputThreshold     = 0,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = F,
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = 1,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE);

  print(paste('done with permutation pass', i))
  PermutatePvalueMatrix[,i]=permuted$cis$min.pv.gene
}

row.names(PermutatePvalueMatrix) <- names(permuted$cis$min.pv.gene)
write.table(t(PermutatePvalueMatrix), permutation_matrix_output_filename, quote=F, sep='\t', col.names = T, row.names=F)



