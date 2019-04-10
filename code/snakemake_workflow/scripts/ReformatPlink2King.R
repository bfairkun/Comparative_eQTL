args = commandArgs(trailingOnly=TRUE)
PlinkKingMatrixIn <- args[1]
KingMatrixOut <- args[2]

nc<-max(count.fields(PlinkKingMatrixIn, sep="\t")) + 1

kingMatrix = rbind(NA, as.matrix(read.table(PlinkKingMatrixIn, sep='\t', col.names=1:nc, fill=T)));
errorCovariance <- Matrix::forceSymmetric(kingMatrix,uplo="L")
errorCovariance[errorCovariance<0] <- 0
errorCovariance[is.na(errorCovariance)] <- 0.5
errorCovariance <- as.matrix(errorCovariance)

write.table(errorCovariance, KingMatrixOut, sep='\t', col.names = F, row.names = F)
