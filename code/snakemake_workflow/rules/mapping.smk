def get_mapped_chunks(wildcards):
    splitid, = glob_wildcards(config["temp_files_prefix"] + "trimmed_chunks/{sample}-{unit}/1.split.{{splitid}}.gz".format(sample=wildcards.sample, unit=wildcards.unit))
    return expand(config["temp_files_prefix"] + "mapped_chunks/{sample}-{unit}/{splitid}.sorted.bam", sample=wildcards.sample, unit=wildcards.unit, splitid=splitid)

rule trim_reads_pe:
    input:
        get_fastq
    output:
        fastq1=temp("trimmed/{sample}-{unit}.1.fastq.gz"),
        fastq2=temp("trimmed/{sample}-{unit}.2.fastq.gz"),
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt"]["pe"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    shell:
        """
        cutadapt {params} -o {output.fastq1} -p {output.fastq2} {input} &> {log}
        """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.split_fastqs_to_chunks.get(**wildcards).output[0]
    return expand(config["temp_files_prefix"] + "mapped_chunks/{sample}-{unit}/{splitid}.sorted.bam",
           sample=wildcards.sample,
           unit=wildcards.unit,
           splitid=glob_wildcards(os.path.join(checkpoint_output, "1.split.{splitid}.gz")).splitid)

# dynamic is being depreciated by snakemake
if config["read_mapping"]["bwa"]:
    rule split_fastqs_to_chunks:
        input:
            # "../../simulatedFastq/1M_Simulated.R1.fastq.gz"
            get_trimmed_reads_tosplit
        output:
            # dynamic("{input[0]}.split.{splitid}")
            R1=dynamic(config["temp_files_prefix"] + "trimmed_chunks/{sample}-{unit}.1.split.{splitid}.gz"),
            R2=dynamic(config["temp_files_prefix"] + "trimmed_chunks/{sample}-{unit}.2.split.{splitid}.gz")
        params:
            split_line_count=config["read_mapping"]["chunksize"]*4
        shell:
            """
            zcat {input[0]} | split -{params.split_line_count} - {config[temp_files_prefix]}trimmed_chunks/{wildcards.sample}-{wildcards.unit}.1.split. && gzip {config[temp_files_prefix]}trimmed_chunks/{wildcards.sample}-{wildcards.unit}.1.split.*
            zcat {input[1]} | split -{params.split_line_count} - {config[temp_files_prefix]}trimmed_chunks/{wildcards.sample}-{wildcards.unit}.2.split. && gzip {config[temp_files_prefix]}trimmed_chunks/{wildcards.sample}-{wildcards.unit}.2.split.*
            """


    rule map_chunks_of_reads_bwa:
        input:
            R1=config["temp_files_prefix"] + "trimmed_chunks/{sample}-{unit}.1.split.{splitid}.gz",
            R2=config["temp_files_prefix"] + "trimmed_chunks/{sample}-{unit}.2.split.{splitid}.gz"
        output:
            bam=temp(config["temp_files_prefix"] + "mapped_chunks/{sample}-{unit}.split.{splitid}.sorted.bam"),
        log:
            "logs/bwa_mem/{sample}-{unit}.split.{splitid}.log"
        params:
            index=config["ref"]["genome"],
            extra=get_read_group
        threads: 8
        shell:
            "(bwa mem -t {threads} {params.extra} {params.index} {input.R1} {input.R2} | samtools sort -o {output.bam}) &> {log}"

    rule merge_bams:
        input:
            dynamic(config["temp_files_prefix"] + "mapped_chunks/{sample}-{unit}.split.{splitid}.sorted.bam")
        output:
            # "trimmed/A-lane1.2.fastq.gz.split.merged"
            temp(config["temp_files_prefix"] + "mapped/{sample}-{unit}.sorted.bam")
        log:
            "logs/samtools_merge/{sample}-{unit}.log"
        shell:
            """
            samtools merge -cf {output} {input} &> {log}
            """

# if config["read_mapping"]["bwa"]:
#     checkpoint split_fastqs_to_chunks:
#         input:
#             get_trimmed_reads_tosplit
#         output:
#             clusters=directory(config["temp_files_prefix"] + "trimmed_chunks/{sample}-{unit}"),
#             # R1=temp(config["temp_files_prefix"] + "trimmed_chunks/{sample}-{unit}/1.split.{splitid}.gz"),
#             # R2=temp(config["temp_files_prefix"] + "trimmed_chunks/{sample}-{unit}/2.split.{splitid}.gz")
#         params:
#             split_line_count=config["read_mapping"]["chunksize"]*4
#         shell:
#             """
#             mkdir -p {config[temp_files_prefix]}trimmed_chunks/{wildcards.sample}-{wildcards.unit} &&
#             zcat {input[0]} | split -{params.split_line_count} - {config[temp_files_prefix]}trimmed_chunks/{wildcards.sample}-{wildcards.unit}/1.split. && gzip {config[temp_files_prefix]}trimmed_chunks/{wildcards.sample}-{wildcards.unit}/1.split.*
#             zcat {input[1]} | split -{params.split_line_count} - {config[temp_files_prefix]}trimmed_chunks/{wildcards.sample}-{wildcards.unit}/2.split. && gzip {config[temp_files_prefix]}trimmed_chunks/{wildcards.sample}-{wildcards.unit}/2.split.*
#             """
#
#
#     rule map_chunks_of_reads_bwa:
#         input:
#             R1=config["temp_files_prefix"] + "trimmed_chunks/{sample}-{unit}/1.split.{splitid}.gz",
#             R2=config["temp_files_prefix"] + "trimmed_chunks/{sample}-{unit}/2.split.{splitid}.gz"
#         output:
#             bam=temp(config["temp_files_prefix"] + "mapped_chunks/{sample}-{unit}/{splitid}.sorted.bam"),
#         log:
#             "logs/bwa_mem/{sample}-{unit}/{splitid}.log"
#         params:
#             index=config["ref"]["genome"],
#             extra=get_read_group
#         threads: 8
#         shell:
#             "(bwa mem -t {threads} {params.extra} {params.index} {input.R1} {input.R2} | samtools sort -o {output.bam}) &> {log}"
#
#
#
#     # rule merge_bams:
#     #     """
#     #     manually removing some temporary files in shell command because snakemake's temp() function does not allow removing directories and the temp() function does not wait for intermediate steps between checkpoint and aggregate rules (because checkpoint output is not listed as input for merge rule)
#     #     """
#     #     input:
#     #         aggregate_input
#     #     output:
#     #         temp(config["temp_files_prefix"] + "mapped/{sample}-{unit}.sorted.bam")
#     #     log:
#     #         "logs/samtools_merge/{sample}-{unit}.log"
#     #     shell:
#     #         """
#     #         samtools merge -cf {output} {input} &> {log} &&
#     #         rm -rf {config[temp_files_prefix]}trimmed_chunks/{wildcards.sample}-{wildcards.unit} &&
#     #         rm -rf {config[temp_files_prefix]}mapped_chunks/{wildcards.sample}-{wildcards.unit}
#     #         """
#
#
#
#     rule merge_bams_from_ancient:
#         """
#         rule is a hack to merge bams when all mapped_chunks are already present, without having to rerun steps with dynamic jobs as snakemake would otherwise do. Mutually exclusive with merge_bams rule. Comment out one or the other
#         """
#         input:
#             #If need to perform mapping as well from completed trimmed chunks, specify get_mapped_chunks as input. If mapping already done, leave input blank.
#             get_mapped_chunks
#         output:
#             temp(config["temp_files_prefix"] + "mapped/{sample}-{unit}.sorted.bam")
#         log:
#             "logs/samtools_merge/{sample}-{unit}.log"
#         shell:
#             """
#             samtools merge -cf {output} {config[temp_files_prefix]}mapped_chunks/{wildcards.sample}-{wildcards.unit}.split.*.sorted.bam &> {log}
#             """

else:
    rule map_reads_hisat:
        input:
            reads=get_trimmed_reads
        output:
            bam=temp(config["temp_files_prefix"] + "mapped/{sample}-{unit}.sorted.bam"),
            qc="qc/hisat2/{sample}-{unit}.txt"
        log:
            "logs/hisat2/{sample}-{unit}.log"
        params:
            index=config["ref"]["genomeindexbase"],
            extra=get_hisat_read_group,
            sort="samtools",
            sort_order="coordinate"
        threads: 8
        shell:
            "(hisat2 -p {threads} {params.extra} --no-spliced-alignment --summary-file {output.qc} --no-softclip --new-summary -x {params.index} -1 {input[0]} -2 {input[1]} | samtools sort -o {output.bam}) &> {log}"


rule mark_duplicates:
    input:
        config["temp_files_prefix"] + "mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=protected("dedup/{sample}-{unit}.sorted.bam"),
        metrics="qc/dedup/{sample}-{unit}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{unit}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    shell:
        "picard MarkDuplicates {params} INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} &> {log}"


rule index_bams:
    input:
        # bam = get_sample_bams_dedup,
        bam = lambda wildcards: expand("dedup/{sample}-{unit}.sorted.bam", sample=wildcards.sample, unit=wildcards.unit)
    output:
        bai="dedup/{sample}-{unit}.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule merge_unit_bams:
    input:
        bam = get_sample_bams_dedup,
    output:
        bam=protected("dedup_merged/{sample}.sorted.bam"),
    shell:
        """
        samtools merge -cf {output.bam} {input.bam}
        """

rule index_merged_bams:
    input:
        bam = lambda wildcards: "dedup_merged/{sample}.sorted.bam".format(sample=wildcards.sample)
    output:
        bai="dedup_merged/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input.bam}"

