Human_genomeDir = config["PowerAnalysis"]["Human_ref"]["genome_dir"][:-1]

rule STAR_make_index_human:
    input:
        fasta=config["PowerAnalysis"]["Human_ref"]["genome_fasta"],
        gtf=config["PowerAnalysis"]["Human_ref"]["genome_gtf"]
    output:
        index = config["PowerAnalysis"]["Human_ref"]["genome_dir"] + "chrLength.txt",
    log:
        "logs/PowerAnalysis/STAR_make_index_human.log"
    params:
        genomeDir = Human_genomeDir
    shell:
        """
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN 4 --genomeDir {params.genomeDir} --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} &> {log}
        """


rule MergeFastqFilesPower:
    input:
        lambda wildcards: PowerAnalysisFastqFrame.at[wildcards.PowerAnalysisSample, "fastq"] 
    output:
        config["temp_files_prefix"] + "PowerAnalysisFastq/{PowerAnalysisSample}.fastq.gz"
    shell:
        "cat {input} > {output}"

rule MergeFastqFilesPowerForTony:
    input:
        lambda wildcards: PowerAnalysisFastqFrame.at[wildcards.PowerAnalysisSample, "fastq"] 
    output:
        "/project2/yangili1/bjf79/FastqForTony/{PowerAnalysisSample}.fastq.gz"
    shell:
        "cat {input} > {output}"

def GetFirstItemOnlyIfNecessary(ListOrString):
    """
    Helper function to reformat outout if pandas df.at method which outputs
    a string if there is only one match but an array of strings if there is
    more than one match... I want the output to always a be string, being the
    first item in the array if necessary.
    """
    if isinstance(ListOrString, str):
        return ListOrString
    else:
        return ListOrString[0]

def GetAlignmentIndexForPowerAnalysis(wildcards):
    Species = GetFirstItemOnlyIfNecessary(PowerAnalysisFastqFrame.at[wildcards.PowerAnalysisSample, "Species"])
    if Species == "Human":
        return Human_genomeDir
    elif Species == "Chimp":
        return "MiscOutput/STAR_index"

def GetSamplesBySpecies(wildcards):
    SampleList = set(PowerAnalysisFastqFrame.index[PowerAnalysisFastqFrame['Species'] == wildcards.species].tolist())
    return expand(config["temp_files_prefix"] + "PowerAnalysisBams/{PowerAnalysisSample}/Aligned.sortedByCoord.out.bam", PowerAnalysisSample = SampleList)


def GetOrthologousExonGtfBySpecies(wildcards):
    if wildcards.species == "Chimp":
        return config["PowerAnalysis"]["Subread"]["Chimp_gtf"]
    elif wildcards.species == "Human":
        return config["PowerAnalysis"]["Subread"]["Human_gtf"]

rule STAR_alignment_PowerAnalysis:
    input:
        index_human = config["PowerAnalysis"]["Human_ref"]["genome_dir"] + "chrLength.txt",
        index_chimp = "MiscOutput/STAR_index/chrLength.txt",
        R1 = config["temp_files_prefix"] + "PowerAnalysisFastq/{PowerAnalysisSample}.fastq.gz"
    log:
        "logs/PowerAnalysis/STAR/{PowerAnalysisSample}.log"
    threads: 8
    params:
        index = GetAlignmentIndexForPowerAnalysis
    output:
        bam = config["temp_files_prefix"] + "PowerAnalysisBams/{PowerAnalysisSample}/Aligned.sortedByCoord.out.bam",
        bai = config["temp_files_prefix"] + "PowerAnalysisBams/{PowerAnalysisSample}/Aligned.sortedByCoord.out.bam.bai",
        star_log = config["temp_files_prefix"] + "PowerAnalysisBams/{PowerAnalysisSample}/Log.final.out"
    shell:
        """
        ulimit -v 31538428657
        STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn <(zcat {input.R1} | fastx_trimmer -l 75) --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate  --outFileNamePrefix {config[temp_files_prefix]}PowerAnalysisBams/{wildcards.PowerAnalysisSample}/ &> {log}
        samtools index {output.bam}
        """

rule PowerAnalysisMultiqc:
    input:
        expand(config["temp_files_prefix"] + "PowerAnalysisBams/{PowerAnalysisSample}/Log.final.out", PowerAnalysisSample=PowerAnalysisFastqFrame.index),
    output:
        "PowerAnalysis/Multiqc/multiqc_report.html"
    shell:
        """
        multiqc -f -o PowerAnalysis/Multiqc/ {input}
        """

rule SubSamplePowerAnalysisBams:
    """
    Before making count table, I want to subsample bams to equal read depth.
    """
    input:
        bam = config["temp_files_prefix"] + "PowerAnalysisBams/{PowerAnalysisSample}/Aligned.sortedByCoord.out.bam",
        bai = config["temp_files_prefix"] + "PowerAnalysisBams/{PowerAnalysisSample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        bam = config["temp_files_prefix"] + "PowerAnalysisBams/{PowerAnalysisSample}/Subsampled/{depth}.bam",
        bai = config["temp_files_prefix"] + "PowerAnalysisBams/{PowerAnalysisSample}/Subsampled/{depth}.bam.bai"
    log:
        "PowerAnalysis/SubSamplePowerAnalysisBams/{PowerAnalysisSample}.{depth}.log"
    shell:
        """
        samtools view -bh -s $(samtools idxstats {input.bam} | awk '$1!="*" {{sum+=$3}} END {{printf "%.4f", {wildcards.depth}/sum}}') {input.bam} > {output.bam} 2> {log}
        samtools index {output.bam}
        """


SedReplace = (config["temp_files_prefix"] + "PowerAnalysisBams/").replace("/", 
"\/")
rule PowerAnalysisSubreadCount:
    input:
        GetSamplesBySpecies
    output:
        Output = "PowerAnalysis/Subread/{species}.subread.txt.gz",
        OutputGit = config["gitinclude_output"] + "PowerAnalysisFullCountTable.{species}.subread.txt.gz",
        OutputSummary = "PowerAnalysis/Subread/{species}.subread.txt.summary"
    log:
        "logs/PowerAnalysis/SubreadCount/{species}.txt"
    params:
        Gtf = GetOrthologousExonGtfBySpecies
    shell:
        """
        ~/software/subread-1.6.3-Linux-x86_64/bin/featureCounts -a {params.Gtf} -o {wildcards.species}.temp.txt {input} &> {log}
        # Format the header for the count table to just contain sample name
        # instead of bam filepath
        cat {wildcards.species}.temp.txt | sed -e '2s/{SedReplace}//g' | sed -e '2s/\/Aligned.sortedByCoord.out.bam//g' | gzip - > {output.Output}
        cat {wildcards.species}.temp.txt.summary | sed -e '1s/{SedReplace}//g' | sed -e '1s/\/Aligned.sortedByCoord.out.bam//g' > {output.OutputSummary}
        rm {wildcards.species}.temp.txt {wildcards.species}.temp.txt.summary
        cp {output.Output} {output.OutputGit}
        """

def GetSamplesBySpeciesDepthMatched(wildcards):
    SampleList = set(PowerAnalysisFastqFrame.index[PowerAnalysisFastqFrame['Species'] == wildcards.species].tolist())
    return expand(config["temp_files_prefix"] + "PowerAnalysisBams/{PowerAnalysisSample}/Subsampled/{depth}.bam", PowerAnalysisSample=SampleList, depth=wildcards.depth)

rule PowerAnalysisSubreadCountDepthMatched:
    input:
        GetSamplesBySpeciesDepthMatched
    output:
        Output = "PowerAnalysis/Subread/{species}.{depth}.subread.txt.gz",
        OutputGit = config["gitinclude_output"] + "PowerAnalysisCountTable.{species}.{depth}.subread.txt.gz",
        OutputSummary = "PowerAnalysis/Subread/{species}.{depth}.subread.txt.summary"
    log:
        "logs/PowerAnalysis/SubreadCount/{species}.{depth}.txt"
    params:
        Gtf = GetOrthologousExonGtfBySpecies
    shell:
        """
        ~/software/subread-1.6.3-Linux-x86_64/bin/featureCounts -a {params.Gtf} -o {wildcards.species}.temp.txt {input} &> {log}
        # Format the header for the count table to just contain sample name
        # instead of bam filepath
        cat {wildcards.species}.temp.txt | sed -e '2s/{SedReplace}//g' | sed -e '2s/\/Subsampled\/{wildcards.depth}.bam//g' | gzip - > {output.Output}
        cat {wildcards.species}.temp.txt.summary | sed -e '1s/{SedReplace}//g' | sed -e '1s/\/Subsampled\/{wildcards.depth}.bam//g' > {output.OutputSummary}
        rm {wildcards.species}.temp.txt {wildcards.species}.temp.txt.summary
        cp {output.Output} {output.OutputGit}
        """

rule DE_edgeR_BootstrapRep:
    input:
        ChimpTable = "../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz",
        HumanTable = "../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz",
        CountTables = GetCountTablesForDE_Bootstrap,
        DropFile = "../../data/DE_SamplesToDrop.txt"
    output:
        Results = "PowerAnalysis/BootstrapReps/{ReadDepthKey}_{seed}.txt.gz"
    log:
        "logs/MakeSupplementalTables/DE_edgeR_BootstrapRep.{ReadDepthKey}.{seed}.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/DE_BootstrapSamples.R {wildcards.seed} {input.CountTables} {wildcards.ReadDepthKey} {output.Results} &> {log}
        """

rule MergePowerBootstrapReps:
    input:
        BootstrapRepOutput
    output:
        "PowerAnalysis/BootstrapRepsMerged.txt.gz"
    shell:
        "cat {input} > {output}"


rule PlotPowerResults:
    input:
        bootsrap = "PowerAnalysis/BootstrapRepsMerged.txt.gz",
        real = "../../output/Final/TableS2.tab",
    output:
        "../../output/Final/FigS_Power_NumDE_{ReadDepthKey}.pdf"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/PlotPowerAnalysisBootstrapResults.R {input.bootsrap} {input.real} {wildcards.ReadDepthKey}
        """


