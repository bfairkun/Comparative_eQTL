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

rule get_tss_and_tata_annotations:
    """
    From paper titled 'refTSS: A Reference Data Set for Human and Mouse
    Transcription Start Sites'
    """
    output:
        tata = "../../data/TSS_annotations/TataAnnotations.txt.gz",
        tss_classifications = "../../data/TSS_annotations/TssClassificationAnnotations.txt.gz",
    shell:
        """
        wget -O {output.tata} http://reftss.clst.riken.jp/datafiles/current/human/tata_box_annotations/hg38_tata_annotation_v3.txt.gz
        wget -O {output.tss_classifications} http://reftss.clst.riken.jp/datafiles/current/human/tss_classification/TSS.classification.hg38.gz
        """

rule SubsetCountTable:
    input:
        CountTable = "OverdispersionAnalysis/FullCountTable.txt.gz",
        SampleList = "../../data/OverdispersionGTExAnalysisSampleList.{list}.txt"
    log:
        "logs/Overdispersion/SubsetCountTable.{list}.log"
    output:
        "OverdispersionAnalysis/TableSubset.{list}.txt"
    shell:
        """
        Rscript scripts/SubsetGTExCountTable.R {input.CountTable} {input.SampleList} {output} &> {log}
        """

rule CalculateOverdispersionTissueMatrix:
    input:
        CountTable = "OverdispersionAnalysis/TableSubset.txt",
        SampleList = "../../data/OverdispersionGTExAnalysisSampleList.{list}.txt",
        GeneLengths="OverdispersionAnalysis/GTEx.genelengths.txt"
    output:
        mu = "OverdispersionAnalysis/GtexTissueMatrix.mu.{list}.txt",
        overdispersion = "OverdispersionAnalysis/GtexTissueMatrix.overdispersion.{list}.txt",
    log:
        "logs/OverdispersionAnalysis/CalculateOverdispersionTissueMatrix.{list}.log"
    shell:
        """
        Rscript scripts/CalculateOverdispersionForGtexTissues.R {input.CountTable} {input.SampleList} {input.GeneLengths} {output.mu} {output.overdispersion} &> {log}
        """

rule CalculateOverdispersionTissueMatrix_NoLengthNorm:
    input:
        CountTable = "OverdispersionAnalysis/TableSubset.txt",
        SampleList = "../../data/OverdispersionGTExAnalysisSampleList.{list}.txt",
        GeneLengths="OverdispersionAnalysis/GTEx.genelengths.txt"
    output:
        mu = "OverdispersionAnalysis/GtexTissueMatrix.mu.NoLengthNorm.{list}.txt",
        overdispersion = "OverdispersionAnalysis/GtexTissueMatrix.overdispersion.NoLengthNorm.{list}.txt",
    log:
        "logs/OverdispersionAnalysis/CalculateOverdispersionTissueMatrix.{list}.log"
    shell:
        """
        Rscript scripts/CalculateOverdispersionForGtexTissues_NoLengthNorm.R {input.CountTable} {input.SampleList} {input.GeneLengths} {output.mu} {output.overdispersion} &> {log}
        """

rule copyToOutput_GTExOverdispersion:
    input:
        mu = "OverdispersionAnalysis/GtexTissueMatrix.mu.NoLengthNorm.{list}.txt",
        overdispersion = "OverdispersionAnalysis/GtexTissueMatrix.overdispersion.NoLengthNorm.{list}.txt",
        mu2 = "OverdispersionAnalysis/GtexTissueMatrix.mu.{list}.txt",
        overdispersion2 = "OverdispersionAnalysis/GtexTissueMatrix.overdispersion.{list}.txt",
    output:
        mu = "../../output/GtexTissueMatrix.mu.NoLengthNorm.{list}.txt.gz",
        overdispersion = "../../output/GtexTissueMatrix.overdispersion.NoLengthNorm.{list}.txt.gz",
        mu2 = "../../output/GtexTissueMatrix.mu.{list}.txt.gz",
        overdispersion2 = "../../output/GtexTissueMatrix.overdispersion.{list}.txt.gz",
    shell:
        """
        cat {input.mu} | gzip - > {output.mu}
        cat {input.overdispersion} | gzip - > {output.overdispersion}
        cat {input.mu2} | gzip - > {output.mu2}
        cat {input.overdispersion2} | gzip - > {output.overdispersion2}
        """

rule OverdispersionInChimps:
    input:
        ChimpCountTable = "../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz",
        HumanCountTable = "../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz",
        Metadata = "../../data/Metadata.xlsx",
        ChimpEgenes = "../../output/ChimpEgenes.eigenMT.txt.gz",
        HumanEgenes = "../../data/GTEX_v8_eGenes/Heart_Left_Ventricle.v8.egenes.txt.gz",
        Biomart = "../../data/Biomart_export.Hsap.Ptro.orthologs.txt.gz"
    output:
        OverdispersionParameters = "../../output/OverdispersionEstimatesFromChimp.txt",
        OverdispersionParameters2 = "../../output/OverdispersionEstimatesFromChimp_NoVirusChallangedIndividuals.txt",
        OverdispersionParameters_NoLengthNorm = "../../output/OverdispersionEstimatesFromChimp.NoLengthNorm.txt",
        OverdispersionParameters2_NoLengthNorm = "../../output/OverdispersionEstimatesFromChimp_NoVirusChallangedIndividuals.NoLengthNorm.txt",
    log:
        "logs/OverdispersionAnalysis/CalculateOverdispersionChimps.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/CalculateChimpOverdispersion.R 2> {log}
        """

rule GetHeartLeftVentricleChromHMM:
    output:
        "../../data/HeartLeftVentricle.ChromHMM.hg38.bed.gz"
    shell:
        "wget -O {output} https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/E095_18_core_K27ac_hg38lift_dense.bed.gz"
