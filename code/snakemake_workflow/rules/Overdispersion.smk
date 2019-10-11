rule get_GTEx_CountTable:
    output:
        "OverdispersionAnalysis/FullCountTable.txt.gz"
    shell:
        """
        wget -O {output} https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
        """

rule get_GTEx_gtf:
    output:
        "OverdispersionAnalysis/GTEx.gencode.v26.gtf"
    shell:
        "wget -O {output} https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf"

rule get_GTEx_geneLength:
    input:
        "OverdispersionAnalysis/GTEx.gencode.v26.gtf"
    output:
        "OverdispersionAnalysis/GTEx.genelengths.txt"
    log:
        "logs/OverdispersionAnalysis/get_GTEx_geneLength.log"
    shell:
        """
        /project2/gilad/bjf79_project1/software/GTFtools_0.6.5/gtftools.py -l {output} {input}
        """

rule SubsetCountTable:
    input:
        CountTable = "OverdispersionAnalysis/FullCountTable.txt.gz",
        SampleList = "../../data/OverdispersionGTExAnalysisSampleList.txt"
    log:
        "logs/Overdispersion/SubsetCountTable.log"
    output:
        "OverdispersionAnalysis/TableSubset.txt"
    shell:
        """
        Rscript scripts/SubsetGTExCountTable.R {input.CountTable} {input.SampleList} {output} &> {log}
        """

rule CalculateOverdispersionTissueMatrix:
    input:
        CountTable = "OverdispersionAnalysis/TableSubset.txt",
        SampleList = "../../data/OverdispersionGTExAnalysisSampleList.txt",
        GeneLengths="OverdispersionAnalysis/GTEx.genelengths.txt"
    output:
        mu = "OverdispersionAnalysis/GtexTissueMatrix.mu.txt",
        overdispersion = "OverdispersionAnalysis/GtexTissueMatrix.overdispersion.txt",
    log:
        "logs/OverdispersionAnalysis/CalculateOverdispersionTissueMatrix"
    shell:
        """
        Rscript scripts/CalculateOverdispersionForGtexTissues.R {input.CountTable} {input.SampleList} {input.GeneLengths} {output.mu} {output.overdispersion} &> {log}
        """