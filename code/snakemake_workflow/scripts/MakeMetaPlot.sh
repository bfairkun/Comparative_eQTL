#!/usr/bin/env bash

set -xe
cat ../eQTL_mapping/MatrixEQTL/Results_lmm/Results.3.txt.top | awk -F'\t' -v OFS='\t' 'NF>1 && $6<2 {split($1,a,":"); print a[1], a[2], a[2] +1, "1"}' | sort | uniq | bedtools sort -i - | awk -F'\t' -v OFS='\t' 'BEGIN {print "track type=bedGraph name=track_label description=center_label"} {print $0}' > ../scratch/test.bg

# cat ../eQTL_mapping/MatrixEQTL/Results_lmm/Results.3.txt | awk -F'\t' 'NF>1 && $6<0.1 {print $2}' | sort | grep -w -f - ../../../data/cDNA.all.chromosomal.bed > ../scratch/test.genes.bed

cat ../eQTL_mapping/MatrixEQTL/Results_lmm/Results.3.txt | awk -F'\t' 'NF>1 {print $2}' | sort | uniq | grep -w -f - ../../../data/cDNA.all.chromosomal.bed > ../scratch/test.genes.bed

bedGraphToBigWig ../scratch/test.bg /project2/gilad/bjf79/genomes/Pan_tro_3.0_Ensembl/Sequence/DNA/Pan_troglodytes.Pan_tro_3.0.dna_sm.chromosome.fa.fai.chromsizes ../scratch/test.bw

computeMatrix scale-regions -o MyOut.gz -R ../scratch/test.genes.bed -S ../scratch/test.bw -b 1000 -a 1000 -bs 50 --unscaled5prime 250 --unscaled3prime 250

plotProfile -m MyOut.gz -out ../scratch/test.meta.png --plotType=fill --averageType sum
