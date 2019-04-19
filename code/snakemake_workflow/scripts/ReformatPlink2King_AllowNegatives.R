args = commandArgs(trailingOnly=TRUE)
PlinkKingMatrixIn <- args[1]
KingMatrixOut <- args[2]

PlinkKingMatrixIn <- "./code/snakemake_workflow/scratch/plink2.king"
PlinkKingMatrixIn <- "./code/snakemake_workflow/scratch/MyKingFilteredGenos.king"
OtherMatrix <- as.matrix(read.table("code/snakemake_workflow/eQTL_mapping/GRM.sXX.txt", sep='\t'))
eigen(OtherMatrix + diag(39)*0.001)$values


nc<-max(count.fields(PlinkKingMatrixIn, sep="\t")) + 1

kingMatrix = rbind(NA, as.matrix(read.table(PlinkKingMatrixIn, sep='\t', col.names=1:nc, fill=T)));
errorCovariance <- Matrix::forceSymmetric(kingMatrix,uplo="L")
errorCovariance[is.na(errorCovariance)] <- 1
errorCovariance <- as.matrix(errorCovariance)
eigen(errorCovariance)$values
svd(errorCovariance)$d

normalize_kinmat <- function(kinmat){
  #normalize kinship so that Kij \in [0,1]
  tmp=kinmat - min(kinmat)
  tmp=tmp/max(tmp)
  tmp[1:9,1:9]
  #fix eigenvalues to positive
  diag(tmp)=diag(tmp)-min(eigen(tmp)$values)
  tmp[1:9,1:9]
  return(tmp)
}
A<-normalize_kinmat(errorCovariance)
dim(A)
eigen(A)$values



write.table(errorCovariance, KingMatrixOut, sep='\t', col.names = F, row.names = F)


mat <- scan("./code/snakemake_workflow/scratch/plink2.king.pythonreformatted.txt")
mat <- matrix(mat, ncol = 39, byrow = TRUE)
