#!/usr/bin/env bash
# Loop through filepaths as commandline arguments

# set -xe #Debug mode

printf "file\tinstrument\trunid\tflowcellid\tflowcelllane\tbarcode"
for MYFILE in "$@"
do
    samplename=$(basename -s ".fastq.gz" $MYFILE)
    printf "\n$MYFILE\t"
    barcode=$(zcat $MYFILE | paste -d'\t'  - - - - | awk -F'\t' '{split($1,a,":"); print a[10]}' | head -10000 | sort | uniq -c | sort -nr | head -1)
    outstr=$(zcat $MYFILE | head -1 | awk -F'\t' -v OFS='\t' '{split($1,a,":"); print a[1], a[2], a[3], a[4]}')
    printf "$outstr\t$barcode"
done
