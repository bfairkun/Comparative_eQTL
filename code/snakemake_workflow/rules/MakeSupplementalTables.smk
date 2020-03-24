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

rule DispersionIterations:
    input:
        ChimpTable = "../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz",
        HumanTable = "../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz",
        DE = "../../output/Final/TableS1.tab",
    output:
        Chimp = "Dispersion/BoostrapSEIterations/Chimp.Chunk.{n}.txt",
        Human = "Dispersion/BoostrapSEIterations/Human.Chunk.{n}.txt",
    params:
        InitialSeed = GetInitialSeedNumberForPermutationChunkBootstrapSE,
        NumberPermutations = config["Overdispersion"]["BootstrapSE_ChunkSize"]
    log:
        "logs/MakeSupplementalTables/DispersionIterationsSE/Chunk.{n}.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/CalculateOverdispersionBoostappedSE.R {input.ChimpTable} {input.HumanTable} {input.DE} {output.Chimp} {output.Human} {params.NumberPermutations} {params.InitialSeed} &> {log}
        """

LastChunkBootstrap_SE=int(config["Overdispersion"]["BootstrapSE_NumChunks"])-1
rule MergeBootstrapChunks:
    input:
        ChimpChunks = expand("Dispersion/BoostrapSEIterations/Chimp.Chunk.{n}.txt", n=range(0, int(config["Overdispersion"]["BootstrapSE_NumChunks"]))),
        HumanChunks = expand("Dispersion/BoostrapSEIterations/Human.Chunk.{n}.txt", n=range(0, int(config["Overdispersion"]["BootstrapSE_NumChunks"]))),
    output:
        Chimp = "Dispersion/BoostrapSE.Chimp.PermutationsCombined.txt",
        Human = "Dispersion/BoostrapSE.Human.PermutationsCombined.txt"
    log:
        "MakeSupplementalTables/MergeBootstrapChunks.log"
    shell:
        """
        cat Dispersion/BoostrapSEIterations/Chimp.Chunk.0.txt > {output.Chimp}
        for i in {{1..{LastChunkBootstrap_SE}}}; do
            tail -n +2 Dispersion/BoostrapSEIterations/Chimp.Chunk.${{i}}.txt >> {output.Chimp}
        done
        cat Dispersion/BoostrapSEIterations/Human.Chunk.0.txt > {output.Human}
        for i in {{1..{LastChunkBootstrap_SE}}}; do
            tail -n +2 Dispersion/BoostrapSEIterations/Human.Chunk.${{i}}.txt >> {output.Human}
        done
        """

rule DispersionIterations_Inference:
    input:
        ChimpTable = "../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz",
        HumanTable = "../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz",
        DE = "../../output/Final/TableS1.tab",
    output:
        Chunks = "Dispersion/BoostrapInferenceIterations/Chunk.{n}.txt",
    params:
        InitialSeed = GetInitialSeedNumberForPermutationChunkBootstrapInference,
        NumberPermutations = config["Overdispersion"]["BoostrapInference_ChunkSize"]
    log:
        "logs/MakeSupplementalTables/DispersionIterationsInference/Chunk.{n}.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/CalculateOverdispersionBoostappedInference.R {input.ChimpTable} {input.HumanTable} {input.DE} {output.Chunks} {params.NumberPermutations} {params.InitialSeed} &> {log}
        """

LastChunkBootstrap_Inference=int(config["Overdispersion"]["BoostrapInference_NumChunks"])-1
rule MergeBootstrapChunks_Inference:
    input:
        Chunks = expand("Dispersion/BoostrapInferenceIterations/Chunk.{n}.txt", n=range(0, int(config["Overdispersion"]["BoostrapInference_NumChunks"]))),
    output:
        Chunks = "Dispersion/BoostrapInference.PermutationsCombined.txt",
    log:
        "MakeSupplementalTables/MergeBootstrapChunks_Inference.log"
    shell:
        """
        cat Dispersion/BoostrapInferenceIterations/Chunk.0.txt > {output.Chunks}
        for i in {{1..{LastChunkBootstrap_Inference}}}; do
            tail -n +2 Dispersion/BoostrapInferenceIterations/Chunk.${{i}}.txt >> {output.Chunks}
        done
        """

rule GetPvaluesFromBoostrapReplicates:
    input:
        Chunks = "Dispersion/BoostrapInference.PermutationsCombined.txt",

rule GetSE_FromBootstrapReplicates:
    input:
        ChimpBootsrapReps = "Dispersion/BoostrapSE.Chimp.PermutationsCombined.txt",
        HumanBootsrapReps = "Dispersion/BoostrapSE.Human.PermutationsCombined.txt"
    output:
        ChimpSEOut = "Dispersion/ChimpSE.tab",
        HumanSEOut = "Dispersion/HumanSE.tab"
    shell:
        """
        ./scripts/GetSE_FromBootstrapIterations.py
        """

rule GetPValues_FromBootstrapReplicates:
    input:
        BootsrapReps = "Dispersion/BoostrapInference.PermutationsCombined.txt",
        Observed = "../../output/OverdispersionEstimatesFromChimp.txt"
    output:
        P = "../../output/OverdispersionEstimatesFromChimp.txt.Pvals.tab"
    shell:
        """
        ./scripts/SubtractBootstrapIterations.py
        """
