
rule DetermineLowcomplexityRegions:
    input:
        ref = config["ref"]["genome"]
    output:
        "MiscOutput/LowComplexityRegions.bed"
    log:
    shell:
        """
        dustmasker -in {input.ref} | awk 'BEGIN {{ OFS='\\t'; sequence="" }} {{ if (match($0, ">.*")) {{ gsub(">", "", $0); gsub(" .*", "", $0); sequence=$0 }} else {{ gsub(" - ", "\\t", $0); print sequence "\\t" $1 "\\t" $2+1 }} }}' > {output}
        """

rule SummarizeFoldCoverageBySample:
    input: expand("qc/samtools-idxstats/{sample}.idxstats", sample=sorted(samples["sample"])),
    output: "MiscOutput/FoldCoveragePerSample.tab"
    shell:
        """
        paste <(awk 'FNR==1 {{print mappedcount, genomelength; mappedcount=0; genomelength=0}} {{mappedcount+=$3; genomelength+=$2}} END {{print mappedcount, genomelength}}' {input}) <(awk 'BEGIN {{print}} FNR==1 {{print FILENAME}}' {input} ) | tail -n +2| awk -v OFS='\\t' '{{print $3, $1*150/$2}}' > {output}
        """

rule DetermineCallableSitesByCoverage_by_contig:
    input:
        bam = expand("dedup_merged/{sample}.sorted.bam", sample=sorted(samples["sample"])),
        bai = expand("dedup_merged/{sample}.sorted.bam.bai", sample=sorted(samples["sample"])),
        FoldCoverageSummary = "MiscOutput/FoldCoveragePerSample.tab" 
    output:
        "MiscOutput/CallableSitesByCoverage/{contig}.bed"
    shell:
        """
        samtools depth -r {wildcards.contig} -d 175 {input.bam} |  python3 scripts/DepthToCallableRegions.py 3 4 {input.FoldCoverageSummary} | bedtools merge -i - > {output}
        """

rule merge_CallableSitesByCoverage_by_contig:
    input:
        bedlist = expand("MiscOutput/CallableSitesByCoverage/{contig}.bed", contig=contigs),
        fai = config["ref"]["genome"] + ".fai"
    output:
        "MiscOutput/CallableSitesByCoverage/FullGenome.bed"
    shell:
        """
        cat {input.bedlist} | bedtools sort -i - -g {input.fai} > {output}
        """

rule DetermineCallableSites:
    input:
        LCRegions = "MiscOutput/LowComplexityRegions.bed",
        AdequateCoverageRegions = "MiscOutput/CallableSitesByCoverage/FullGenome.bed",
        fai = config["ref"]["genome"] + ".fai"
    output:
        "MiscOutput/CallableSites.bed"
    shell:
        """
        bedtools subtract -sorted -g {input.fai} -a {input.AdequateCoverageRegions} -b {input.LCRegions} > {output}
        """

rule CountCallableSites:
    input:
        Callable = "MiscOutput/CallableSites.bed",
        fai = config["ref"]["genome"] + ".fai"
    output:
        "MiscOutput/CallableSites.summary"
    shell:
        """
        printf "callable sites: $(awk -F'\\t' '{{sum+=$3-$2}} END{{print sum}}' {input.Callable})\\n" > {output}
        printf "total sites: $(awk -F'\\t' '{{sum+=$2}} END{{print sum}}' {input.fai})\\n" >> {output}
        """

rule FilterReadsAtUncallableSites:
    """
    Some regions contain extremely high coverage at regions that are uncallable for variant calling. These regions are filtered from bam to reduce memory footprint for freebayes. More specifically, reads are left in if they are in CallableRegions (which excludes these high coverage regions)
    """
    input:
        CallableSites = "MiscOutput/CallableSites.bed",
        bam = lambda wildcards: "dedup_merged/{sample}.sorted.bam".format(sample=wildcards.sample),
        bai = lambda wildcards: "dedup_merged/{sample}.sorted.bam.bai".format(sample=wildcards.sample),
    output:
        bam = config["temp_files_prefix"] + "dedupAndCallableBams/{sample}.sorted.bam"
    log:
        "logs/bedtools/IntersectingCallable_{sample}.log"
    shell:
        """
        bedtools intersect -a {input.bam} -b {input.CallableSites} -u -sorted > {output.bam}
        """
rule index_filtered_bams:
    input:
        bam = lambda wildcards: config["temp_files_prefix"] + "dedupAndCallableBams/{sample}.sorted.bam".format(sample=wildcards.sample)
    output:
        bai=config["temp_files_prefix"] + "dedupAndCallableBams/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"


if config["variant_calling"]["freebayes_in_chunks"]:
    rule call_variants_freebayes_chunks:
        input:
            bam = expand([config["temp_files_prefix"] + "dedupAndCallableBams/{sample}.sorted.bam"], sample=samples["sample"]),
            bai = expand([config["temp_files_prefix"] + "dedupAndCallableBams/{sample}.sorted.bam.bai"], sample=samples["sample"]),
            ref=config["ref"]["genome"],
        output:
            vcfgz="genotyped/all.{freebayes_region_chunk}.vcf"
        log:
            "logs/freebayes/{freebayes_region_chunk}.log"
        params:
            extra="--min-coverage 3 --max-coverage 150 -k --standard-filters --report-genotype-likelihood-max -n 2",
            intervals=get_freebayes_chunked_region
        shell:
            "freebayes {params.intervals} {params.extra} -f {input.ref} -b {input.bam} > {output} 2> {log}"

    rule merge_freebayes_chunks_by_contig:
        input:
            vcf=get_vcf_chunks_by_contig
        output:
            vcfgz="genotyped/ByContig/all.{contig}.vcf.gz"
        log:
            "logs/freebayes/merging.{contig}.vcf.log"
        params:
            SortMem="20G"
        shell:
            "(bcftools concat {input} | bcftools sort -m {params.SortMem} -O z > {output}) &> {log}"

else:
    rule call_variants_freebayes:
        input:
            bam = expand(["dedupAndCallableBams/{u.sample}-{u.unit}.sorted.bam"], u=units.itertuples()),
            bai = expand(["dedupAndCallableBams/{u.sample}-{u.unit}.sorted.bam.bai"], u=units.itertuples()),
            ref=config["ref"]["genome"],
            known=config["ref"]["known-variants"]
        output:
            vcfgz="genotyped/ByContig/all.{contig}.vcf.gz"
        log:
            "logs/freebayes/{contig}.log"
        params:
            extra="--prob-contamination 0.05 -k --no-population-priors",
            intervals="-r {contig}"
        shell:
            "(freebayes {params.intervals} {params.extra} -f {input.ref} -b {input.bam} > {output}) 2> {log}"

rule IndexInitCalls:
    input: lambda wildcards: "genotyped/ByContig/all.{contig}.vcf.gz".format(contig=wildcards.contig)
    output: "genotyped/ByContig/all.{contig}.vcf.gz.tbi"
    shell: "tabix -p vcf {input}"

rule FilterInitialCalls:
    """need to add bcftools norm filter because through some bug (https://github.com/samtools/bcftools/issues/221) bcftools view with -R will output some duplicate lines. It's also good practice to normalize calls anywy"""
    input:
        InitCalls=ancient (lambda wildcards: "genotyped/ByContig/all.{contig}.vcf.gz".format(contig=wildcards.contig) ),
        CallableSites= ancient( "MiscOutput/CallableSites.bed"),
        Tabix= ancient (lambda wildcards: "genotyped/ByContig/all.{contig}.vcf.gz.tbi".format(contig=wildcards.contig) ),
    output:
        Filtered="filtered/{contig}.vcf.gz",
        Tabix="filtered/{contig}.vcf.gz.tbi"
    log: "logs/FilteringInitCalls/{contig}.log"
    params:
        VariantFilter = "-i '%QUAL>30'",
        SortMem = "20G",
        fasta = config["ref"]["genome"]
    shell:
        """
        (bcftools view {params.VariantFilter} -R {input.CallableSites}  {input.InitCalls} | bcftools norm -f {params.fasta} -d none | bcftools sort -m {params.SortMem} -O z > {output.Filtered}) 2> {log}
        tabix -p vcf {output.Filtered} &>> {log}
        """

rule MergeFilteredCalls:
    input:
        FilteredCalls=expand("filtered/{contig}.vcf.gz", contig=contigs),
        Tabix=expand("filtered/{contig}.vcf.gz.tbi", contig=contigs),
    output:
        vcf = "filtered/all.vcf.gz",
        tbi = "filtered/all.vcf.gz.tbi"
    log:
        "logs/FilteringInitCalls/Merge.log"
    params:
        SortMem="20G"
    shell:
        """
        (bcftools concat {input.FilteredCalls} | bcftools sort -O z > {output.vcf} ) &> {log}
        tabix -p vcf {output.vcf} &>> {log}
        """

# rule IndexFilteredCalls:
#     input: "filtered/all.vcf.gz"
#     output: "filtered/all.vcf.gz.tbi"
#     log: "logs/FilteringInitCalls/Merge.index.log"
#     shell:
#         "tabix -p vcf {input} &> {log}"

rule bcftools_stats:
    input:
        "filtered/all.vcf.gz"
    output:
        "MiscOutput/filtered.bcftools.stats"
    log:
        "logs/Misc/bcftools_stats.log"
    shell:
        """
        bcftools stats {input} > {output} 2> {log}
        """

rule bcftools_smplstats:
    input:
        "filtered/all.vcf.gz"
    output:
        "MiscOutput/filtered.bcftools.smplstats"
    log:
        "logs/Misc/bcftools_smplstats"
    shell:
        """
        bcftools +smpl-stats {input} > {output} 2> {log}
        """

