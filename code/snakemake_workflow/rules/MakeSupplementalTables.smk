rule DE_edgeR:
    input:
        ChimpTable = "../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz",
        HumanTable = "../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz",
        GeneIDs = "../../data/HumanGeneIdToHGNC_Symbol.Biomart.txt.gz"
    output:
        "../../output/Final/TableS1.tab"
    log:
        "logs/MakeSupplementalTables/DE_edgeR.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/DE_edgeR.R &> {log}
        """
