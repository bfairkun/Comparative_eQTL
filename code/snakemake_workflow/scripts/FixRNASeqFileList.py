#!/usr/bin/env python3

DNASeqFileListSamples = set()
with open("samples.tsv", 'rU') as f:
    f.readline()
    for line in f:
        sample = line.strip('\n')
        DNASeqFileListSamples.add(sample)
DNASeqDict = {sample.upper():sample for sample in DNASeqFileListSamples}
print DNASeqDict

RNASeqFileListSamples = set()
with open("RNASeqFileList.tsv", 'rU') as f:
    with open("RNASeqFileList.corrected.tsv", 'w') as fout:
        fout.write("sample\tfastq\n")
        f.readline()
        for line in f:
            sample, fastq = line.strip('\n').split('\t')
            RNASeqFileListSamples.add(sample)
            fout.write("{newsample}\t{fastq}\n".format(fastq=fastq, newsample=DNASeqDict[sample.upper()]))
