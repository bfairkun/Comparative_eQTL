# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)
library(tidyverse)
# library(qvalue)

# args = commandArgs(trailingOnly=TRUE)
# SNP_file_name <- args[1]
# snps_location_file_name <- args[2]
# expression_file_name <- args[3]
# gene_location_file_name <- args[4]
# covariates_file_name <- args[5]
# errorCovariance_file <- args[6]
# output_file_name_cis <- args[7]
# ouput_QQ <- args[8]
# cisDistance <- args[9]

# setwd("/project2/gilad/bjf79_project1/projects/Comparative_eQTL/code/snakemake_workflow/")
# SNP_file_name <- "eQTL_mapping/SharedPolymorphisms/SpeciesSharedSnps.PanTro5.snps"
# snps_location_file_name <- "eQTL_mapping/SharedPolymorphisms/SpeciesSharedSnps.PanTro5.snploc"
# expression_file_name <- "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt"
# gene_location_file_name <- "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt"
# covariates_file_name <- "../../output/Covariates/0GenotypePCs_and_10RNASeqPCs.covariates"
# errorCovariance_file <- "eQTL_mapping/Kinship/GRM.cXX.txt"
# output_file_name_cis = tempfile()
# cisDistance = 1E6
# SNP_file_name_matched <- "eQTL_mapping/SharedPolymorphisms/SpeciesSharedMatchedSnps.PanTro5.snps"
# snps_location_file_name_matched <- "eQTL_mapping/SharedPolymorphisms/SpeciesSharedMatchedSnps.PanTro5.snploc"
# OutCisImage <- "~/temp/Cis.png"
# OutTransImage <- "~/temp/Trans.png"

setwd("/project2/gilad/bjf79_project1/projects/Comparative_eQTL/code/snakemake_workflow/")
SNP_file_name <- "GTEX_renalysis/LeflerSnps/MatrixEQTL/Genotypes.txt"
snps_location_file_name <- "GTEX_renalysis/LeflerSnps/MatrixEQTL/snploc.txt"
expression_file_name <- "GTEX_renalysis/LeflerSnps/MatrixEQTL/expression.txt"
gene_location_file_name <- "GTEX_renalysis/LeflerSnps/MatrixEQTL/geneloc.txt"
covariates_file_name <- "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Heart_Left_Ventricle.v8.covariates.txt"
errorCovariance_file <- "eQTL_mapping/Kinship/GRM.cXX.txt"
output_file_name_cis = tempfile()
cisDistance = 1E6
SNP_file_name_matched <- "GTEX_renalysis/LeflerSnps/MatrixEQTL/Genotypes.control.txt"
snps_location_file_name_matched <- "GTEX_renalysis/LeflerSnps/MatrixEQTL/snploc.control.txt"
OutCisImage <- "../../output/LeflerTestedSnps.Human.Cis.QQ.png"
OutTransImage <- "../../output/LeflerTestedSnps.Human.Trans.QQ.png"
CisResults <- "../../output/LeflerTestedSnps.Human.cis.tsv"

# SNP_file_name <- args[1]
# snps_location_file_name <- args[2]
# expression_file_name <- args[3]
# gene_location_file_name <- args[4]
# covariates_file_name <- args[5]
# errorCovariance_file <- args[6]
# output_file_name_cis = tempfile()
# cisDistance = args[7]
# SNP_file_name_matched <- args[8]
# snps_location_file_name_matched <- args[9]
# OutCisImage <- args[10]
# OutTransImage <- args[11]

# SNP_file_name <- "code/snakemake_workflow/eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps"
# snps_location_file_name <- "code/snakemake_workflow/eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc"
# expression_file_name <- "code/snakemake_workflow/eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt"
# gene_location_file_name <- "code/snakemake_workflow/eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt"
# covariates_file_name <- "output/Covariates/0GenotypePCs_and_10RNASeqPCs.covariates"
# errorCovariance_file <- "code/snakemake_workflow/eQTL_mapping/Kinship/GRM.cXX.txt"
# output_file_name_cis = tempfile()
# cisDistance = 250000

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS


output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance <- as.matrix(read.table(errorCovariance_file,sep='\t'))


# Distance for local gene-SNP pairs
cisDist = as.numeric(cisDistance);
print(cisDist)

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

snps_matched = SlicedData$new();
snps_matched$fileDelimiter = "\t";      # the TAB character
snps_matched$fileOmitCharacters = "NA"; # denote missing values;
snps_matched$fileSkipRows = 1;          # one row of column labels
snps_matched$fileSkipColumns = 1;       # one column of row labels
snps_matched$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps_matched$LoadFile(SNP_file_name_matched);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}



## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
snpspos_matched = read.table(snps_location_file_name_matched, header = TRUE, stringsAsFactors = FALSE);

genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);



#### Real pass
#Include all pvalues in output so that qvalues can be calculated. Filter after.
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = NULL,
  pvOutputThreshold     = 1,
  useModel = useModel,
  verbose = TRUE,
  output_file_name.cis = NULL,
  pvOutputThreshold.cis = 1,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

#Pass with matched snps
me_matched = Matrix_eQTL_main(
  snps = snps_matched,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = NULL,
  pvOutputThreshold     = 1,
  useModel = useModel,
  verbose = TRUE,
  output_file_name.cis = NULL,
  pvOutputThreshold.cis = 1,
  snpspos = snpspos_matched,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
print('done with real pass')

#### make data for permutation pass
ActualData.ExpressionMatrix <- read.table(expression_file_name, header=T, row.names = 1, check.names = FALSE )
ActualData.Cov <- read.table(covariates_file_name, header=T, row.names = 1, check.names = FALSE )

### Run permutation pass
print("Running permutation pass ")

# Permute column labels for both (using same seed for randomization)
# technically I am permuting the column data, and preserving the column labels. That way, I do not have to do the same for the (huge) genotype file, which would take more computational time
set.seed(0)
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
  pvOutputThreshold     = 1,
  useModel = useModel,
  verbose = F,
  output_file_name.cis = NULL,
  pvOutputThreshold.cis = 1,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

print('done with permutation pass')

plot(me)
plot(me_matched)
plot(permuted)

length(me$trans$eqtls$pvalue)

CisQQPlot <- ggplot(me$cis$eqtls, aes(y=-log10(sort(pvalue)), x=-log10(1:length(pvalue)/length(pvalue)))) +
  geom_point(aes(color="Variants shared with human")) +
  geom_point(data=permuted$cis$eqtls, aes(color="Variants shared with human; permuted data")) +
  geom_point(data=me_matched$cis$eqtls, aes(color="Matched control variants")) +
  xlab("-log10(Theoretical-Pvalues)") +
  ylab("-log10(Observed-Pvalues)") +
  geom_abline() +
  scale_color_manual(values = c("Variants shared with human" = "red",
                                "Variants shared with human; permuted data" = "black",
                                "Matched control variants" = "blue") ) +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  theme(legend.key.width=unit(0.2, "cm"),
        legend.direction = "vertical")

ggsave(OutCisImage, plot=CisQQPlot, height=3.7, width=3.7)

wilcox.test(me_matched$cis$eqtls$pvalue, me$cis$eqtls$pvalue)
wilcox.test(me_matched$trans$eqtls$pvalue, me$trans$eqtls$pvalue)


TransQQPlot <- ggplot(me$trans$eqtls, aes(y=-log10(sort(pvalue)), x=-log10(1:length(pvalue)/length(pvalue)))) +
  geom_point(aes(color="Variants shared with chimp")) +
  geom_point(data=permuted$trans$eqtls, aes(color="Variants shared with chimp; permuted data")) +
  geom_point(data=me_matched$trans$eqtls, aes(color="Matched control variants")) +
  xlab("-log10(Theoretical-Pvalues)") +
  ylab("-log10(Observed-Pvalues)") +
  geom_abline() +
  scale_color_manual(values = c("Variants shared with chimp" = "red",
                                "Variants shared with chimp; permuted data" = "black",
                                "Matched control variants" = "blue") ) +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  theme(legend.key.width=unit(0.2, "cm"),
        legend.direction = "vertical")

ggsave(OutTransImage, plot=TransQQPlot, height=3.7, width=3.7)

write.table(me$cis$eqtls, file = CisResults, quote=F, row.names = F, sep='\t')
