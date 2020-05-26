rule CopyNovelRNASeqData:
    input:
        lambda wildcards: GEORNASeqBasenameToFastq[wildcards.basename]
    output:
        "DataForSubmission/GEO/Fastq/{basename}"
    shell:
        """
        cp {input} {output}
        """

rule CheckSumRNASeq:
    input:
        RNASeqFastqList = expand( "DataForSubmission/GEO/Fastq/{basename}", basename=GEORNASeqBasenameToFastq.keys())
    output:
        "DataForSubmission/GEO/FastqChecksums.txt"
    shell:
        """
        md5sum {input} > {output}
        """

rule MetaDataFormatHelper:
    input:
        "DataForSubmission/GEO/FastqChecksums.txt",
        "../../data/NovelRNASeqFiles.tsv"
    output:
        "DataForSubmission/GEO/MetadataHelper.tsv"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript 
        """

rule CopyProcessedFiles:
    input:
        ChimpTable = "PowerAnalysis/Subread/Chimp.subread.txt.gz",
        HumanTable = "PowerAnalysis/Subread/Human.subread.txt.gz",
        ChimpRawTable = "RNASeq/STAR/CountTable.txt.gz",
        ChimpNormalizedTable = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt"
    output:
        ChimpTable = "DataForSubmission/GEO/OrthologousExonRawGeneCounts_Chimpanzee.tsv.gz",
        HumanTable = "DataForSubmission/GEO/OrthologousExonRawGeneCounts_Human.tsv.gz",
        ChimpRawTable = "DataForSubmission/GEO/ChimpanzeesRawGeneCount_eQTL_Mapping.tsv.gz",
        ChimpNormalizedTable = "DataForSubmission/GEO/ChimpanzeesQQNorm_eQTL_Mapping.tsv.gz"
    shell:
        """
        cp {input.ChimpTable} {output.ChimpTable}
        cp {input.HumanTable} {output.HumanTable}
        cp {input.ChimpRawTable} {output.ChimpRawTable}
        cat {input.ChimpNormalizedTable} | gzip - > {output.ChimpNormalizedTable}
        """
