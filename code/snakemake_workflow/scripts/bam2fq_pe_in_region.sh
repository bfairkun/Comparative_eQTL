#!/usr/bin/env bash
if [ "$1" == "-h" ]; then
    echo "Usage: `basename $0`  <chr:Start-Stop> <in.bam> <R1_filename.fastq> <R2_filename.fastq> <Other samtools bam2fq parameters>"
    echo ""
    echo "similar to samtools bam2fq, to outputs bam alignments to paired end read fastq files, but with required parameter of region. Will output both pairs of reads if and only if both pairs of reads map within given region, to produce R1 and R2 files of equal length"
    exit 0
fi

set -xe #Debug mode

region=$1
bam=$2
R1=$3
R2=$4
shift 4
bam2fqflags=$*

cat <(samtools view -H $bam) <(samtools view -f2 -F2304 $bam $region | grep -f <(samtools view -f2 -F2304 $bam $region | awk -F'\t' '{print $1}' | sort | uniq -c | awk '$1==2 {print $2}')) | samtools bam2fq - -1 $R1 -2 $R2
