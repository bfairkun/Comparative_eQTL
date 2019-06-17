import pandas as pd
import os
from snakemake.utils import validate
from collections import defaultdict
from itertools import chain

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
samples_ForMockData = pd.read_table(config["MockDataMaker"]["SampleListForMockData"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

# contigs in reference genome
contigs = pd.read_table(config["ref"]["genome"] + ".fai",
                        header=None, usecols=[0], squeeze=True, dtype=str)

if config["variant_calling"]["freebayes_in_chunks"]:
    freebayes_region_chunks = pd.read_table(os.popen("bedtools makewindows -w " + str(config["variant_calling"]["chunksize"]) + " -g " + config["ref"]["genome"] + ".fai | awk -F'\t' '{print $1\".\"$2\"-\"$3}' "), header=None, squeeze=True, dtype=str, usecols=[0])
    freebayes_region_chunks_dict = defaultdict(list)
    for entry in freebayes_region_chunks:
        freebayes_region_chunks_dict[entry.split(".")[0]].append(entry.split(".")[1])
else:
    freebayes_region_chunks = ""

RNASeqFastqList = pd.read_table(config["RNASeqFileList"], squeeze=True, dtype=str).set_index(["sample"])
RNASeqSampleToFastq_dict = defaultdict(list)
RNASeqBasenameToFastq = dict()
with open(config["RNASeqFileList"]) as RNASeqFileList_fh:
    RNASeqFileList_fh.readline()
    for line in RNASeqFileList_fh:
        samplename, filepath = line.strip('\n').split('\t')
        RNASeqSampleToFastq_dict[samplename].append(filepath)
        RNASeqBasenameToFastq[os.path.basename(filepath)] = filepath

PowerAnalysisFastqFrame = pd.read_csv(config["PowerAnalysis"]["RNASeqFileList"],sep='\t', index_col=0, comment='#')

##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),
    contig="|".join(contigs),
    freebayes_region_chunk="|".join(freebayes_region_chunks),
    freebayes_chunk_coords="|".join([i.split(".")[1] for i in freebayes_region_chunks]),
    species="Chimp|Human"


##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    return + units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def getFastqFromBase(wildcards):
    return(RNASeqBasenameToFastq[wildcards.fastq_basename])

def get_RNASeq_merged_fastq(wildcards):
    """Get merged RNA fastq file for a given {sample} wildcard"""
    return config["temp_files_prefix"] + "RNASeqFastq/{sample}.fastq.gz".format(sample=wildcards.RNASeqSample)

def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def is_single_end2(sample, unit,splitid):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{prefix}{sample}\tSM:{prefix}{sample}\tPL:{platform}'".format(
        prefix=config["read_mapping"]["ReadGroupID_prefix"],
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])

def get_hisat_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return "--rg ID:{sample} --rg SM:{prefix}{sample} --rg PL:{platform}".format(
        prefix=config["read_mapping"]["ReadGroupID_prefix"],
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

def get_trimmed_reads_tosplit(wildcards):
    """Get trimmed reads of given sample-unit."""
    return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                     group=[1, 2], **wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)

def get_sample_bams_dedup(wildcards):
    """Get all aligned reads of given sample."""
    return expand("dedup/{sample}-{unit}.sorted.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)

def get_freebayes_chunked_region(wildcards):
    """Return a string to mark files by the region used in freeybayes command"""
    return "-r " + wildcards.freebayes_region_chunk.replace('.',':')


def get_vcf_chunks_by_contig(wildcards):
    """get_vcf_chunks_by_contig"""
    return expand("genotyped/all.{{contig}}.{coords}.vcf",
                  contig=wildcards.contig,
                  coords=freebayes_region_chunks_dict[wildcards.contig])

genotypedregions, = glob_wildcards("genotyped/all.{region}.vcf.gz")

