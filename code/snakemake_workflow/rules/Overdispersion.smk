rule get_GTEx_CountTable:
    output:
        "OverdispersionAnalysis/FullCountTable.txt.gz"
    shell:
        """
        wget -O {output} https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
        """

rule SubsetCountTable:
    input: "OverdispersionAnalysis/FullCountTable.txt.gz"
    log: "logs/"
    output: "OverdispersionAnalysis/TableSubset.txt.gz"
