#!/usr/bin/env bash

set -xe

# bg for all eqtls
cat  ../eQTL_mapping/MatrixEQTL/Results_lmm_QQNorm_250kB/Results.0GenotypePCs_and_10RNASeqPCs.covariates.txt | awk -F'\t' -v OFS='\t' 'NF>1 && $6<0.3 {split($1,a,"."); print a[2], a[3], a[3] +1, "1"}' | sort | uniq | bedtools sort -i - | awk -F'\t' -v OFS='\t' 'BEGIN {print "track type=bedGraph name=track_label description=center_label"} {print $0}' > ../scratch/test.bg

# bg for best snp-gene pairs
# cat ../eQTL_mapping/MatrixEQTL/Results_lmm_QQNorm_250kB/Results.0GenotypePCs_and_10RNASeqPCs.covariates.txt | awk -F'\t' -v OFS='\t' 'NF>1 && $6<0.3 {print $2,$6,$1}' | sort -k1,1n -k2,2g | awk '!a[$1] {a[$1] =$2} $2 ==a[$1]' | awk -v OFS='\t' '{split($3,a,"."); print a[2], a[3], a[3] +1, "1"}' | sort | uniq | bedtools sort -i - | awk -F'\t' -v OFS='\t' 'BEGIN {print "track type=bedGraph name=track_label description=center_label"} {print $0}' > ../scratch/test.bg

# # bg for all tested snps
# awk -F'\t' -v OFS='\t' 'NR>1 && $2!=23 {print $2,$3,$3+1,0.001}' ../eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc | sort -k1,1 -k2,2n | awk -v OFS='\t' 'BEGIN {print "track type=bedGraph name=track_label description=center_label"} {print $0}' > ../scratch/test.snps.bg


# cat  ../eQTL_mapping/MatrixEQTL/Results_lmm_QQNorm_250kB/Results.0GenotypePCs_and_10RNASeqPCs.covariates.txt | awk -F'\t' 'NF>1 && $6<0.3 {print $2}' | sort | uniq | grep -w -f - ../../../output/Genes.bed > ../scratch/test.genes.bed


# bedGraphToBigWig ../scratch/test.bg /project2/gilad/bjf79/genomes/Pan_tro_3.0_Ensembl/Sequence/DNA/Pan_troglodytes.Pan_tro_3.0.dna_sm.chromosome.fa.fai.chromsizes ../scratch/test.bw
# bedGraphToBigWig ../scratch/test.snps.bg /project2/gilad/bjf79/genomes/Pan_tro_3.0_Ensembl/Sequence/DNA/Pan_troglodytes.Pan_tro_3.0.dna_sm.chromosome.fa.fai.chromsizes ../scratch/test.snps.bw

computeMatrix scale-regions -o MyOut.gz -R ../scratch/test.genes.bed  -S ../scratch/test.bw ../scratch/test.snps.bw -b 250000 -m 50000 -a 250000 -bs 2500

# computeMatrix reference-point -o MyOut.gz -R ../scratch/test.genes.bed -S ../scratch/test.bw -b 500000 -a 500000 -bs 50000 --referencePoint TES

plotProfile -m MyOut.gz -out ../scratch/test.meta.png --plotType=fill --averageType sum
