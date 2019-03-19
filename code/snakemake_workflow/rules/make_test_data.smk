# Given already aligned fastq files for all samples, subsample bam in specified
# region and output fastq files which can serve assd create useful test data.
# 
# How to use: If snakemake is executed with CollectAllMockData as target rule
# (snakemake CollectAllMockData), then the below rules will be applied. Then,
# based on that sample list and the region defined in the config file, mock
# fastq files will be created and a accompanying units.tsv and
# RNASeqFileList.tsv.  I use these files to test the snakemake pipeline using
# smaller datasets, usually by executing the snakemake with the --directory
# parameter in a serperate directory after moving the appropriate units.tsv and
# RNAseqSampleList.tsv files (generated here) to that directory.


RegionStringReformatted = config["MockDataMaker"]["Region"].replace(":",".").replace(",","")

rule CollectAllMockData:
    input:
        expand("test_data/WGS/{sample}." + RegionStringReformatted + ".1.fastq.gz", sample=samples_ForMockData["sample"]),
        # expand("test_data/RNASeq/{sample}." + RegionStringReformatted + ".1.fastq.gz", sample=samples_ForMockData["sample"]),

rule MakeWGS_Fastq:
    input:
        bam = lambda wildcards:"dedup_merged/{sample}.sorted.bam".format(sample=wildcards.sample),
        bai = lambda wildcards:"dedup_merged/{sample}.sorted.bam.bai".format(sample=wildcards.sample),
    output:
        R1 = "test_data/WGS/{sample}." + RegionStringReformatted + ".1.fastq.gz",
        R2 = "test_data/WGS/{sample}." + RegionStringReformatted + ".2.fastq.gz",
    log:
        "logs/makedata/WGS/{sample}.log"
    params:
        region= config["MockDataMaker"]["Region"]
    shell:
        """
        bash scripts/bam2fq_pe_in_region.sh {params.region} {input.bam} {output.R1} {output.R2} &> {log}
        """

rule MakeRNASeq_Fastq:
    input:
        bam = lambda wildcards:"RNASeq/STAR/{sample}/Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample),
        bai = lambda wildcards:"RNASeq/STAR/{sample}/Aligned.sortedByCoord.out.bam.bai".format(sample=wildcards.sample)
    output:
        R1 = "test_data/RNASeq/{sample}." + RegionStringReformatted + ".1.fastq.gz",
        R2 = "test_data/RNASeq/{sample}." + RegionStringReformatted + ".2.fastq.gz",
    log:
        "logs/makedata/RNAseq/{sample}.log"
    params:
        region=config["MockDataMaker"]["Region"]
    shell:
        """
        bash scripts/bam2fq_pe_in_region.sh {params.region} {input.bam} {output.R1} {output.R2} &> {log}
        """
# 
# Not finished
# rule MakeTsvFiles:
#     input:
#     output:
#         unitsTsv = "test_data/units.tsv",
#         RNASeqTsv = "test_data/RNASeqSampleList.tsv"
#     log:
#         ""
#     params:
#         ""
#     run:
#         """
#         with open(output.unitsTsv, "w") as out:
#             for
#         with open(output.RNASeqTsv, "w") as out:
#             for sample in 
#         """


