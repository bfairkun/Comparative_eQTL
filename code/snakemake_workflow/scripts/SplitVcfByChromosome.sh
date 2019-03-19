#!/usr/bin/env sh

set -xe
VcfGzIn=$1
Prefix=$2

# cd $(mktemp -d)
pushd $(mktemp -d)
tempdir=$(pwd)

zcat $VcfGzIn | awk '{chrom=substr($1,1,1); if (chrom=="#"){chrom="header"}; print > chrom ".vcf" }'

shopt -s extglob
for filename in !(header.vcf):
do
    echo $filename
    cat header.vcf $filename > "reheader_$filename"
    bgzip "reheader_$filename"
done
shopt -u extglob

popd
for filename in $tempdir:
do
    echo $filename
done
