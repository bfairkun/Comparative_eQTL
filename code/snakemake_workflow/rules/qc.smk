# rule fastqc_RNA_seq:
#     input:
#         lambda wildcards: config["temp_files_prefix"] + "RNASeqFastq/{RNASeqSample}.fastq.gz".format(RNASeqSample=wildcards.RNASeqSample)
#     output:
#         html="qc/fastqc/RNASeq/{RNASeqSample}.html",
#         zip="qc/fastqc/RNASeq/{RNASeqSample}.zip"
#     shell:
#     wrapper:
#         "0.27.1/bio/fastqc"


rule fastqc:
    input:
        get_fastq
    output:
        html="qc/fastqc/{sample}-{unit}.html",
        zip="qc/fastqc/{sample}-{unit}.zip"
    # shell:
    wrapper:
        "0.27.1/bio/fastqc"

rule samtools_idxstats:
    input:
        lambda wildcards: "dedup_merged/{sample}.sorted.bam".format(sample=wildcards.sample)
    output:
        "qc/samtools-idxstats/{sample}.idxstats"
    shell:
        "samtools idxstats {input} > {output}"

rule qualimap:
    input:
        get_sample_bams_dedup
    output:
        "qc/qualimap/{sample}-{unit}/qualimapReport.html"
    log:
        "logs/qualimap/{sample}-{unit}.log"
    shell:
        """
        qualimap bamqc -bam {input}  -nt 12 -outdir qc/qualimap/{wildcards.sample}-{wildcards.unit} --java-mem-size=29G &> {log}
        """

rule QC_to_publish:
    input:
        bcftools_smplstats = "MiscOutput/filtered.bcftools.smplstats",
        Coverage = "MiscOutput/FoldCoveragePerSample.tab",
        RNASeqCov = "scratch/multiqcOnRNASeq/multiqc_data/multiqc_general_stats.txt"
    output:
        bcftools_smplstats = "../../output/QC/bcftools.smplstats.tab",
        Coverage = "../../output/QC/FoldCoverPerSample.tab",
        RNASeqCov = "../../output/QC/RNASeqMultiQC.stats.tab"
    shell:
        """
        cp {input.bcftools_smplstats} {output.bcftools_smplstats}
        cp {input.Coverage} {output.Coverage}
        cp {input.RNASeqCov} {output.RNASeqCov}
        """

#rule multiqc:
#    input:
#        expand(["qc/qualimap/{u.sample}-{u.unit}/qualimapReport.html",
#                "qc/fastqc/{u.sample}-{u.unit}.zip",
#                "qc/samtools-idxstats/{u.sample}.idxstats",
#                "qc/dedup/{u.sample}-{u.unit}.metrics.txt"],
#               u=units.itertuples()),
#        expand(["RNASeq/STAR/{RNASeqSample}/ReadsPerGene.out.tab",
#                "RNASeq/STAR/{RNASeqSample}/Log.final.out",
#                "qc/fastqc/RNASeq/{RNASeqSample}.zip"],
#                RNASeqSample=RNASeqSampleToFastq_dict.keys()),
#        #"snpeff/all.csv"
#    output:
#        report("qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
#    log:
#        "logs/multiqc.log"
#    shell:
#    wrapper:
#        "0.27.1/bio/multiqc"
