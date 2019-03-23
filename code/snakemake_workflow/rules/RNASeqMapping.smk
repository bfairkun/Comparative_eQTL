
rule MergeFastqFiles:
    input:
        lambda wildcards: expand("{RNASeqSample}", RNASeqSample=RNASeqSampleToFastq_dict[wildcards.RNASeqSample])
    output:
        config["temp_files_prefix"] + "RNASeqFastq/{RNASeqSample}.fastq.gz"
    shell:
        "cat {input} > {output}"

rule STAR_make_index:
    input:
        fasta=config["ref"]["genome"],
        gtf=config["ref"]["genomegtf"]
    output:
        "MiscOutput/STAR_index/chrLength.txt",
    log:
        "logs/STAR/MakingIndex.log"
    shell:
        """
        mkdir -p "MiscOutput/STAR_index/"
        STAR --runMode genomeGenerate --runThreadN 4 --genomeDir MiscOutput/STAR_index/ --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} &> {log}
        """

rule STAR_alignment:
    input:
        index = "MiscOutput/STAR_index/chrLength.txt",
        fastq = get_RNASeq_merged_fastq
    log:
        "logs/STAR/{RNASeqSample}.log"
    threads: 12
    output:
        "RNASeq/STAR/{RNASeqSample}/ReadsPerGene.out.tab",
        "RNASeq/STAR/{RNASeqSample}/Aligned.sortedByCoord.out.bam",
        "RNASeq/STAR/{RNASeqSample}/Log.final.out"
    shell:
        "STAR --runThreadN {threads} --genomeDir MiscOutput/STAR_index/ --readFilesIn {input.fastq} --outSAMtype BAM SortedByCoordinate --outWigStrand Unstranded --outWigType wiggle --alignEndsType EndToEnd --quantMode GeneCounts --twopassMode Basic --readFilesCommand zcat --outFileNamePrefix RNASeq/STAR/{wildcards.RNASeqSample}/ &> {log}"

rule index_RNA_seq_bams:
    input:
        bam = lambda wildcards: "RNASeq/STAR/{sample}/Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample)
    output:
        bai= "RNASeq/STAR/{sample}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"

rule kallisto_make_index:
    input:
        config["ref"]["transcriptsfa"]
    output:
        "MiscOutput/kallisto.idx"
    log:
        "logs/kallisto/MakeIndex.log"
    shell:
        "kallisto index -i {output} {input} &> {log}"

rule kallisto_quant:
    input:
        index = "MiscOutput/kallisto.idx",
        fastq = get_RNASeq_merged_fastq
    output:
        "RNASeq/kallisto/{RNASeqSample}/abundance.tsv"
    log:
        "logs/kallisto/{RNASeqSample}.log",
    shell:
        "kallisto quant -i {input.index} -o RNASeq/kallisto/{wildcards.RNASeqSample}/ --single -l 180 -s 50 {input.fastq} &> {log}"

rule kallisto_quant_merge_into_count_table:
    input:
        quant=expand("RNASeq/kallisto/{sample}/abundance.tsv", sample=RNASeqSampleToFastq_dict.keys()),
    output:
        config["gitinclude_output"] + "CountTable.tpm.txt"
    shell:
        """
        paste <(awk '{{print $1}}' {input.quant[0]}) \
        <(cat <(echo $(ls -1v {input.quant}) | sed 's/RNASeq\/kallisto\///g' | sed 's/\/abundance.tsv//g' | sed 's/[[:blank:]]+/\t/g') <(awk '{{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 }} END {{ for(i=1;i<=FNR;i++) print a[i] }}' $(ls -1v {input.quant}) | awk -v OFS='\t' 'NR>1')) > {output}
        """

rule kallisto_quant_each_fastq:
    input:
        index = "MiscOutput/kallisto.idx",
        # fastq = lambda(wildcards): RNASeqBasenameToFastq[wildcards.fastq_basename]
        fastq = getFastqFromBase
    output:
        "RNASeq/kallisto_per_fastq/{fastq_basename}/abundance.tsv"
    log:
        "logs/kallisto_per_fastq/{fastq_basename}.log"
    shell:
        "kallisto quant -i {input.index} -o RNASeq/kallisto_per_fastq/{wildcards.fastq_basename}/ --single -l 180 -s 50 {input.fastq} &> {log}"

rule kallisto_quant_merge_into_count_table_each_fastq:
    input:
        quant=expand("RNASeq/kallisto_per_fastq/{fastq_basename}/abundance.tsv", fastq_basename=RNASeqBasenameToFastq.keys())
    output:
        config["gitinclude_output"] + "CountTable.SeparatedByFastq.tpm.txt"
    shell:
        """
        paste <(awk '{{print $1}}' {input.quant[0]}) \
        <(cat <(echo $(ls -1v {input.quant}) | sed 's/RNASeq\/kallisto\///g' | sed 's/\/abundance.tsv//g' | sed 's/[[:blank:]]+/\t/g') <(awk '{{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 }} END {{ for(i=1;i<=FNR;i++) print a[i] }}' $(ls -1v {input.quant}) | awk -v OFS='\t' 'NR>1')) > {output}
        """
rule Gather_RNA_seq_FlowCellInfo:
    input: expand("{myfiles}", myfiles=RNASeqBasenameToFastq.values())
    output: "../../output/RNA_seq_FlowCellInfo.txt"
    log: "logs/Gather_RNA_seq_FlowCellInfo"
    shell: "bash scripts/GetFastqIdentifierInfo.sh {input}"
# rule feature_counts:
#     input:
#         bam=expand("RNASeq/STAR/{sample}/Aligned.sortedByCoord.out.bam", sample=RNASeqSampleToFastq_dict.keys())
#     output:
#         "CountTable.txt"
#     shell:
#         "echo {input.bam}"
