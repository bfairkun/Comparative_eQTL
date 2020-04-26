rule DE_edgeR:
    input:
        ChimpTable = "../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz",
        HumanTable = "../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz",
        GeneIDs = "../../data/HumanGeneIdToHGNC_Symbol.Biomart.txt.gz"
    output:
        "../../output/Final/TableS2.tab"
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
        DE = "../../output/Final/TableS2.tab",
    output:
        Chimp = "Dispersion/BoostrapSEIterations/Chimp{DropSampleList}.Chunk.{n}.txt",
        Human = "Dispersion/BoostrapSEIterations/Human{DropSampleList}.Chunk.{n}.txt",
    params:
        InitialSeed = GetInitialSeedNumberForPermutationChunkBootstrapSE,
        NumberPermutations = config["Overdispersion"]["BootstrapSE_ChunkSize"],
        DropSampleFile = lambda wildcards: list(samples_dispersion[samples_dispersion['WildcardIdentifier']==wildcards.DropSampleList]['DropSampleList'])[0]
    wildcard_constraints:
        DropSampleList = "|".join(list(samples_dispersion["WildcardIdentifier"]))
    log:
        "logs/MakeSupplementalTables/DispersionIterationsSE/Chunk.{n}.{DropSampleList}.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/CalculateOverdispersionBoostappedSE.R {input.ChimpTable} {input.HumanTable} {input.DE} {output.Chimp} {output.Human} {params.NumberPermutations} {params.InitialSeed} {params.DropSampleFile} &> {log}
        """


LastChunkBootstrap_SE=int(config["Overdispersion"]["BootstrapSE_NumChunks"])-1
rule MergeBootstrapChunks:
    input:
        ChimpChunks = expand("Dispersion/BoostrapSEIterations/Chimp{{DropSampleList}}.Chunk.{n}.txt", n=range(0, int(config["Overdispersion"]["BootstrapSE_NumChunks"]))),
        HumanChunks = expand("Dispersion/BoostrapSEIterations/Human{{DropSampleList}}.Chunk.{n}.txt", n=range(0, int(config["Overdispersion"]["BootstrapSE_NumChunks"]))),
    output:
        Chimp = "Dispersion/BoostrapSE.Chimp{DropSampleList}.PermutationsCombined.txt",
        Human = "Dispersion/BoostrapSE.Human{DropSampleList}.PermutationsCombined.txt"
    log:
        "MakeSupplementalTables/MergeBootstrapChunks.{DropSampleList}.log"
    wildcard_constraints:
        DropSampleList = "|".join(list(samples_dispersion["WildcardIdentifier"]))
    shell:
        """
        cat Dispersion/BoostrapSEIterations/Chimp{wildcards.DropSampleList}.Chunk.0.txt > {output.Chimp}
        for i in {{1..{LastChunkBootstrap_SE}}}; do
            tail -n +2 Dispersion/BoostrapSEIterations/Chimp{wildcards.DropSampleList}.Chunk.${{i}}.txt >> {output.Chimp}
        done
        cat Dispersion/BoostrapSEIterations/Human{wildcards.DropSampleList}.Chunk.0.txt > {output.Human}
        for i in {{1..{LastChunkBootstrap_SE}}}; do
            tail -n +2 Dispersion/BoostrapSEIterations/Human{wildcards.DropSampleList}.Chunk.${{i}}.txt >> {output.Human}
        done
        """

rule DispersionIterations_Inference:
    input:
        ChimpTable = "../../output/PowerAnalysisFullCountTable.Chimp.subread.txt.gz",
        HumanTable = "../../output/PowerAnalysisFullCountTable.Human.subread.txt.gz",
        DE = "../../output/Final/TableS2.tab",
    output:
        Chunks = "Dispersion/BoostrapInferenceIterations/Chunk{DropSampleList}.{n}.txt",
    params:
        InitialSeed = GetInitialSeedNumberForPermutationChunkBootstrapInference,
        NumberPermutations = config["Overdispersion"]["BoostrapInference_ChunkSize"],
        DropSampleFile = lambda wildcards: list(samples_dispersion[samples_dispersion['WildcardIdentifier']==wildcards.DropSampleList]['DropSampleList'])[0]
    wildcard_constraints:
        DropSampleList = "|".join(list(samples_dispersion["WildcardIdentifier"]))
    log:
        "logs/MakeSupplementalTables/DispersionIterationsInference/Chunk{DropSampleList}.{n}.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/CalculateOverdispersionBoostappedInference.R {input.ChimpTable} {input.HumanTable} {input.DE} {output.Chunks} {params.NumberPermutations} {params.InitialSeed} {params.DropSampleFile} &> {log}
        """

LastChunkBootstrap_Inference=int(config["Overdispersion"]["BoostrapInference_NumChunks"])-1
rule MergeBootstrapChunks_Inference:
    input:
        Chunks = expand("Dispersion/BoostrapInferenceIterations/Chunk{{DropSampleList}}.{n}.txt", n=range(0, int(config["Overdispersion"]["BoostrapInference_NumChunks"]))),
    output:
        Chunks = "Dispersion/BoostrapInference.{DropSampleList}PermutationsCombined.txt",
    wildcard_constraints:
        DropSampleList = "|".join(list(samples_dispersion["WildcardIdentifier"]))
    log:
        "MakeSupplementalTables/MergeBootstrapChunks_Inference.{DropSampleList}.log"
    shell:
        """
        cat Dispersion/BoostrapInferenceIterations/Chunk{wildcards.DropSampleList}.0.txt > {output.Chunks}
        for i in {{1..{LastChunkBootstrap_Inference}}}; do
            tail -n +2 Dispersion/BoostrapInferenceIterations/Chunk{wildcards.DropSampleList}.${{i}}.txt >> {output.Chunks}
        done
        """


rule GetSE_FromBootstrapReplicates:
    input:
        ChimpBootsrapReps = "Dispersion/BoostrapSE.Chimp{DropSampleList}.PermutationsCombined.txt",
        HumanBootsrapReps = "Dispersion/BoostrapSE.Human{DropSampleList}.PermutationsCombined.txt"
    output:
        ChimpSEOut = "Dispersion/Chimp{DropSampleList}SE.tab",
        HumanSEOut = "Dispersion/Human{DropSampleList}SE.tab"
    wildcard_constraints:
        DropSampleList = "|".join(list(samples_dispersion["WildcardIdentifier"]))
    shell:
        """
        ./scripts/GetSE_FromBootstrapIterations.py {input.ChimpBootsrapReps} {input.HumanBootsrapReps} {output.ChimpSEOut} {output.HumanSEOut}
        """

rule CopySE_toOutput:
    input:
        ChimpSEOut = "Dispersion/Chimp{DropSampleList}SE.tab",
        HumanSEOut = "Dispersion/Human{DropSampleList}SE.tab"
    output:
        "../../output/OverdispersionEstimatesFromChimp{DropSampleList}.txt.SE.tab.gz"
    wildcard_constraints:
        DropSampleList = "|".join(list(samples_dispersion["WildcardIdentifier"]))
    shell:
        """
        paste -d'\\t' {input.ChimpSEOut} {input.HumanSEOut} | awk -F'\\t' -v OFS='\\t' 'NR==1 {{ print "gene", "Chimp.SE", "Human.SE" }} NR>1 {{ print $1, $2, $4 }}' | gzip - > {output}
        """

rule GetPValues_FromBootstrapReplicates:
    input:
        BootsrapReps = "Dispersion/BoostrapInference.PermutationsCombined.txt",
        Observed = lambda wildcards: list(samples_dispersion[samples_dispersion['WildcardIdentifier']==wildcards.DropSampleList]['ObservedDispersionOutput'])[0]
    output:
        P = "../../output/OverdispersionEstimatesFromChimp{DropSampleList}.txt.Pvals.tab"
    wildcard_constraints:
        DropSampleList = "|".join(list(samples_dispersion["WildcardIdentifier"]))
    shell:
        """
        ./scripts/SubtractBootstrapIterations.py {input.BootsrapReps} {input.Observed} {output.P}
        """

rule AddQvals_CombinedTable:
    input:
        genes = "../../data/HumanGeneIdToHGNC_Symbol.Biomart.txt.gz",
        SE ="../../output/OverdispersionEstimatesFromChimp{DropSampleList}.txt.SE.tab.gz",
        Observed = lambda wildcards: list(samples_dispersion[samples_dispersion['WildcardIdentifier']==wildcards.DropSampleList]['ObservedDispersionOutput'])[0],
        P = "../../output/OverdispersionEstimatesFromChimp{DropSampleList}.txt.Pvals.tab"
    wildcard_constraints:
        DropSampleList = "|".join(list(samples_dispersion["WildcardIdentifier"]))
    output:
        "../../output/Overdispersion{DropSampleList}_P_SE_Combined.txt.gz"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/AddQValue.R {input.Observed} {input.SE} {input.P} {output}
        """

