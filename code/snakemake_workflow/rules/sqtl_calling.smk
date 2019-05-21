# include: "common.smk"
# rule all:
#     input:
#         expand ("sQTL_mapping/juncfiles/{sample}.junc", sample=samples["sample"] ),
#         "sQTL_mapping/leafcutter/juncfilelist.txt",
#         "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz",

rule STAR_to_leafcutter_junc:
    input:
        "RNASeq/STAR/{sample}/SJ.out.tab"
    output:
        "sQTL_mapping/juncfiles/{sample}.junc"
    log:
        "logs/sQTL_mapping/STAR_to_leafcutter_junc/{sample}.log"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$4==1 && $1!="MT" {{ print $1,$2,$3,".",$7,"+" }} $4==2&& $1!="MT" {{ print $1,$2,$3,".",$7,"-" }}' {input} > {output}
        """

rule make_leafcutter_juncfile:
    input:
        expand ("sQTL_mapping/juncfiles/{sample}.junc", sample=samples["sample"] ),
    output:
        "sQTL_mapping/leafcutter/juncfilelist.txt"
    params:
        SamplesToRemove = "MD_And"
    run:
        import os
        with open(output[0], "w") as out: 
            for filepath in input:
                samplename = os.path.basename(filepath).split(".junc")[0]
                if samplename != params.SamplesToRemove:
                    out.write(filepath + '\n')

rule leafcutter_cluster:
    input:
        "sQTL_mapping/leafcutter/juncfilelist.txt",
    output:
        "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz",
        "sQTL_mapping/leafcutter/clustering/leafcutter_perind_numers.counts.gz"
    log:
        "logs/sQTL_mapping/leafcutter_cluster.log"
    shell:
        """
        leafcutter_cluster.py -j {input} -r sQTL_mapping/leafcutter/clustering/ &> {log}
        """

rule leafcutter_prepare_phenotype_table:
    input:
        counts = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz",
        blacklist_chromosomes = "scratch/chromsomseblacklist.txt"
    output:
        phenotypes = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_Catted.txt",
        PCs = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.PCs"
    shell:
        """
        ~/miniconda3/bin/python2.7 ~/CurrentProjects/leafcutter/scripts/prepare_phenotype_table.py -p 13 --ChromosomeBlackList {input.blacklist_chromosomes}  {input.counts}

        cat <(head -1 sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_chr3) <(awk 'FNR>1' sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_chr*) > {output.phenotypes}
        """

rule prepare_MatrixEQTL_for_sQTL:
    input:
        phenotypes = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_Catted.txt",
        fam = config['gitinclude_output'] + "ForAssociationTesting.temp.fam",
    output:
        phenotypes = "sQTL_mapping/MatrixEQTL/sQTL_phenoTable.txt",
        phenotypesReordered = "sQTL_mapping/MatrixEQTL/sQTL_phenoTable.Reordered.txt",
        intron_locs = "sQTL_mapping/MatrixEQTL/sQTL_intron.locs"
    shell:
        """
        cut -d $'\\t' -f 4- {input.phenotypes} > {output.phenotypes}
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "gene", "chr", "start", "stop" }} NR>1 {{ split($4,a,":"); print $4, a[1], a[2], a[3] }}' {input.phenotypes} > {output.intron_locs}
        Rscript scripts/ReorderPhenotypeTableForMatrixEQTL.R {output.phenotypes} {input.fam} {output.phenotypesReordered}
        """

MatrixSQTLModels = expand("sQTL_mapping/MatrixEQTL/Results/Results.PCs.{B}.covariates.txt", B= list(range (1, config["eQTL_mapping"]["CovariatePCs"]["RNASeqPC_max"] + 1) ))
MatrixSQTLCovariates = expand("sQTL_mapping/Covariates/FromLeafcutter.PCs.{M}.covariates.txt", M=list(range (1, config["eQTL_mapping"]["CovariatePCs"]["RNASeqPC_max"] + 1)))

rule make_covariate_file_sQTL:
    """
    Use R to make covariate file that matches order of fam file
    """
    input:
        fam = config['gitinclude_output'] + "ForAssociationTesting.temp.fam",
        leafcutterPCs = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.PCs",
    output:
        MatrixSQTLCovariates
    shell:
        """
        Rscript scripts/SelectNtoM_PCs_toCovariateFiles.R  {input.leafcutterPCs} {input.fam} sQTL_mapping/Covariates/FromLeafcutter.PCs. 0 {config[eQTL_mapping][CovariatePCs][RNASeqPC_max]}
        """

rule MatrixEQTL_sQTL:
    """Matrix EQTL script performs one cis-eqtl scan with real data and one
    scan with permutated sample labels for phenotypes for an empirical null."""
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        phenotypes = "sQTL_mapping/MatrixEQTL/sQTL_phenoTable.Reordered.txt",
        gene_loc = "sQTL_mapping/MatrixEQTL/sQTL_intron.locs",
        covariates = "sQTL_mapping/Covariates/FromLeafcutter.PCs.{covariate_set}.covariates.txt",
        GRM = CovarianceMatrix,
    output:
        results = "sQTL_mapping/MatrixEQTL/Results/Results.PCs.{covariate_set}.covariates.txt",
        fig = "sQTL_mapping/MatrixEQTL/Results/images/Results.{covariate_set}.png",
        permutated_fig = "sQTL_mapping/MatrixEQTL/Results/images/PermutatedResults.{covariate_set}.png",
        permuted_results = "sQTL_mapping/MatrixEQTL/Results/PermutatedResults.{covariate_set}.txt",
    log:
        "logs/sQTL_mapping/MatrixEQTL/{covariate_set}.log"
    shell:
        """
        Rscript scripts/MatrixEqtl_Cis.R {input.snps} {input.snp_locs} {input.phenotypes} {input.gene_loc} {input.covariates} {input.GRM} {output.results} {output.fig} {output.permuted_results} {output.permutated_fig} 100000 &> {log}
        """

rule PlotPCsVsSQTLs:
    input:
        MatrixSQTLModels
    output:
        CattedResult = "sQTL_mapping/MatrixEQTL/Results/ConcatenatedResult.txt",
        Plot = "sQTL_mapping/MatrixEQTL/Results/images/PCsVsSQTLs.pdf"
    log:
        "logs/sQTL_mapping/PlotPCsVsEQTLs.log"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' 'FNR>1 && $6<0.3 {{ print $1,$2,$3,$4,$5,$6,FILENAME  }}' {input} > {output.CattedResult}
        Rscript scripts/Plot_EQTLs_vs_PCs.R {output.CattedResult} {output.Plot}
        """


