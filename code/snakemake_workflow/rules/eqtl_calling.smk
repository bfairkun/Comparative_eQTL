# rule all:
#     input:
#         "../../output/Genes.bed",

# Use covariate files in folder specified in config, otherwise, use all the
# covariates definied by the min and max number of PCs specified in the config
# Note functions must be defined above the rules that use them
if config["eQTL_mapping"]["CovariatesDir"]:
    def get_covariates(wildcards):
        """Get covariate file"""
        return config["eQTL_mapping"]["CovariatesDir"] + "{covariate_set}".format(covariate_set = wildcards.covariate_set)
    CovariateFiles, = glob_wildcards(config["eQTL_mapping"]["CovariatesDir"] + "{CovariateName}")
    MatrixEQTLModels = expand("eQTL_mapping/MatrixEQTL/Results/Results.{CovariateSetName}.txt", CovariateSetName = CovariateFiles)
else:
    def get_covariates(wildcards):
        return "../../output/Covariates/{covariate_set}".format(covariate_set = wildcards.covariate_set)
    MatrixEQTLModels = expand("eQTL_mapping/MatrixEQTL/Results/Results.{A}GenotypePCs_and_{B}RNASeqPCs.covariates.txt", 
        A= list( range( config["eQTL_mapping"]["CovariatePCs"]["GenotypePC_min"], config["eQTL_mapping"]["CovariatePCs"]["GenotypePC_max"] + 1) ),
        B= list( range( config["eQTL_mapping"]["CovariatePCs"]["RNASeqPC_min"], config["eQTL_mapping"]["CovariatePCs"]["RNASeqPC_max"] + 1) ),
    )

rule make_plink_file_from_vcf_for_testing:
    input:
        vcf = ancient("PopulationSubstructure/ReferencePanelMerged.annotated.vcf.gz"),
    output:
        "eQTL_mapping/plink/Unfiltered.bed"
    log:
        "logs/eQTL_mapping/make_plink_file_from_vcf_for_testing.log"
    shell:
        """
        plink --id-delim '-' --vcf {input.vcf} --vcf-half-call m --allow-extra-chr --make-bed --out eQTL_mapping/plink/Unfiltered &> {log}
        """

rule filter_plink_file_for_testing:
    """
    .fam file output to the gitinclude output directory for convenience, so
    that I can easily access it from my local laptop with git, to edit scripts
    that choose covariates and such with my local RStudio... For the covariate
    files must have samples ordered same as the fam (phenotype) file
    """
    input:
        bed = "eQTL_mapping/plink/Unfiltered.bed"
    output:
        bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        fam = config['gitinclude_output'] + "ForAssociationTesting.temp.fam",
        bim = "eQTL_mapping/plink/ForAssociationTesting.bim"
    log:
        "logs/eQTL_mapping/filter_plink_file_for_testing.log"
    params:
        stdparams = "--keep-fam eQTL_mapping/plink/KeepFam.txt --memory 28000",
        extra = '--remove eQTL_mapping/plink/Remove.txt --maf 0.1 --geno --hwe 1e-7.5'
    shell:
        """
        plink --bfile eQTL_mapping/plink/Unfiltered  --allow-extra-chr --make-bed --out eQTL_mapping/plink/ForAssociationTesting.temp {params.stdparams} {params.extra} &> {log}
        mv eQTL_mapping/plink/ForAssociationTesting.temp.bed {output.bed}
        mv eQTL_mapping/plink/ForAssociationTesting.temp.bim {output.bim}
        mv eQTL_mapping/plink/ForAssociationTesting.temp.fam {output.fam}
        """

rule make_plink_file_for_GRM:
    """Note that GRM will ignore individuals with a missing phenotype -9 in the fam file"""
    input:
        bed = "eQTL_mapping/plink/Unfiltered.bed"
    output:
        bed = "eQTL_mapping/Kinship/ForGRM.bed",
        fam = "eQTL_mapping/Kinship/ForGRM.fam"
    params:
        stdparams = "--keep-fam eQTL_mapping/plink/KeepFam.txt",
        extra = '--remove eQTL_mapping/plink/Remove.txt --maf 0.05 --geno --hwe 1e-7.5'
    shell:
        """
        plink --bfile eQTL_mapping/plink/Unfiltered  --allow-extra-chr --make-bed --out eQTL_mapping/Kinship/ForGRM {params.stdparams} {params.extra}
        sed -i 's/-9$/1/' {output.fam}
        """

rule GetGeneLevelLocation:
    input: config["ref"]["genomegtf"]
    output: "../../output/Genes.bed"
    shell:
        """
        cat {input} | awk -F '\\t' -v OFS='\\t' '$3=="gene"' | perl -lne '/^(.+?)\\t.+?\\t.+?\\t(.+?)\\t(.+?)\\t.+?\\t(.+?)\\t.+?\\tgene_id "(.+?)".+$/ and print "$1\\t$2\\t$3\\t$5\\t.\\t$4"' > {output}
        """


rule make_cis_gene_windows:
    input:
        genes_bed = "../../output/Genes.bed"
    output:
        "eQTL_mapping/Misc/cis_gene_windows.bed"
    log:
        "logs/eQTL_mapping/cis_gene_windows.log"
    params:
        MaxDistFromTSS = config["eQTL_mapping"]["TSS_flanking_snp_width"],
        faidx = config["ref"]["genome"] + ".fai",
    shell:
        """
        (bedtools sort -i {input.genes_bed} | bedtools flank -l 1 -r 0 -s -g {params.faidx} | bedtools slop -b {params.MaxDistFromTSS} -g {params.faidx} | bedtools sort -i - > {output}) &> {log}
        """

# Different rule for if filtering phenotype matrix if using kallisto vs STAR. For example, with STAR I require 80% of samples to have >10 reads in a gene. In kallisto this is not possible so I have different critera for filtering.
if config["eQTL_mapping"]["read_mapper"] == "kallisto":

    rule kallisto_counttable_to_log10TPM_matrix:
        """
        Use R to filter TPM count matrix for genes to map eqtls for. Additionally,
        do quantile normalization here.
        """
        input:
            fam = config['gitinclude_output'] + "ForAssociationTesting.temp.fam",
            counttable = "RNASeq/kallisto/CountTable.tpm.txt.gz",
            genes_bed = "../../output/Genes.bed",
            transcripts_to_genes = "../../data/Biomart_export.Pan_Tro_3.geneids.txt"
        output:
            log10TPM = config["gitinclude_output"] + "ExpressionMatrix.un-normalized.txt.gz"
        log:
            "logs/eQTL_mapping/kallisto_counttable_to_log10TPM_matrix.log"
        shell:
            """
            Rscript ../../analysis/20190327_MakeFamPhenotypeFile.R {input.counttable} {input.fam} {config[gitinclude_output]}ExpressionMatrix.un-normalized.txt  {input.genes_bed} {input.transcripts_to_genes} &> {log}
            gzip {config[gitinclude_output]}ExpressionMatrix.un-normalized.txt
            """

elif config["eQTL_mapping"]["read_mapper"] == "STAR":

    rule STAR_counttable_to_log10CPM_matrix:
        """
        Use R to filter TPM count matrix for genes to map eqtls for. Additionally,
        do quantile normalization here.
        """
        input:
            fam = config['gitinclude_output'] + "ForAssociationTesting.temp.fam",
            counttable = "RNASeq/STAR/CountTable.txt.gz",
            genes_bed = "../../output/Genes.bed",
        output:
            log10TPM = config["gitinclude_output"] + "ExpressionMatrix.un-normalized.txt.gz"
        log:
            "logs/eQTL_mapping/STAR_counttable_to_log10CPM_matrix.log"
        shell:
            """
            Rscript scripts/STARRawCountTableToLog10CPM.R {input.counttable} {input.fam} {config[gitinclude_output]}ExpressionMatrix.un-normalized.txt {input.genes_bed} &> {log}
            gzip {config[gitinclude_output]}ExpressionMatrix.un-normalized.txt
            """

rule quantile_normalize:
    input: config["gitinclude_output"] + "ExpressionMatrix.un-normalized.txt.gz"
    output: config["gitinclude_output"] + "ExpressionMatrix.normalized.txt.gz"
    log: "logs/eQTL_mapping/quantile_normalize.log"
    shell:
        """
        python3 scripts/StandardizeAndQuantileNormalize.py {input} {output} &>> {log}
        """

if config["eQTL_mapping"]["quantile_normalize"]:
    ExpressionMatrix = config["gitinclude_output"] + "ExpressionMatrix.normalized.txt.gz"
else:
    ExpressionMatrix = config["gitinclude_output"] + "ExpressionMatrix.un-normalized.txt.gz"
rule convert_expression_matrix_to_fam:
    input:
        ExpressionMatrix = ExpressionMatrix,
        fam = config['gitinclude_output'] + "ForAssociationTesting.temp.fam",
    output:
        gene_list = "eQTL_mapping/plink/ForAssociationTesting.fam.phenotype-list",
        fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
    log:
        "logs/eQTL_mapping/convert_expression_matrix_to_fam.log"
    shell:
        """
        Rscript scripts/MergeFamWithPhenotypes.R {input.ExpressionMatrix} {input.fam} {output.fam} {output.gene_list} &>> {log}
        """

rule prune_plink_files_for_GRM:
    input:
        bed = "eQTL_mapping/Kinship/ForGRM.bed",
    output:
        bed =  "eQTL_mapping/Kinship/ForAssociationTesting.pruned.bed",
        fam = "eQTL_mapping/Kinship/ForAssociationTesting.pruned.fam"
    log:
        "logs/eQTL_mapping/prune_plink_files.log"
    shell:
        """
        plink --bfile eQTL_mapping/Kinship/ForGRM --allow-extra-chr --indep-pairwise 50 5 0.5 &> {log}
        plink --bfile eQTL_mapping/Kinship/ForGRM  --allow-extra-chr --extract plink.prune.in --make-bed --out eQTL_mapping/Kinship/ForAssociationTesting.pruned &> {log}
        sed -i 's/-9$/1/' {output.fam}
        rm plink.prune.in plink.prune.out
        """

rule Make_GRM:
    """Genetetic relatedness matrix for gemma LMM"""
    input:
        bed =  "eQTL_mapping/Kinship/ForAssociationTesting.pruned.bed",
    output:
        GRM = "eQTL_mapping/Kinship/GRM.cXX.txt",
        GRM_ouput = "../../output/GRM.cXX.txt"
    log:
        "logs/eQTL_mapping/make_GRM.log"
    shell:
        """
        gemma -gk 1 -bfile eQTL_mapping/Kinship/ForAssociationTesting.pruned -o GRM -outdir eQTL_mapping/Kinship &> {log}
        cp {output.GRM} {output.GRM_output}
        """

rule Make_King_GRM:
    """Genetic relatedness matrix from plink2 (king robust)"""
    input:
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
    output:
        GRM = "eQTL_mapping/Kinship/GRM.king.txt"
    log:
        "logs/eQTL_mapping/make_King_GRM.log"
    shell:
        """
        plink2 --bfile eQTL_mapping/plink/ForAssociationTesting --king
        """

rule Subset_GenesToTest:
    input:
        genes_to_test = "eQTL_mapping/plink/ForAssociationTesting.fam.phenotype-list",
    output:
        subset_to_test = "eQTL_mapping/GeneListToTest.txt",
    params:
        n = 20000
    shell:
        """
        awk -v OFS='\\t' '{{ print NR, $1 }}' {input.genes_to_test} | shuf -n {params.n} > {output.subset_to_test}
        """

rule GetCisSNPsJobList:
    input:
        genes_bed = "eQTL_mapping/Misc/cis_gene_windows.bed",
        subset_to_test = "eQTL_mapping/GeneListToTest.txt",
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
    output:
        batchshellscript = "eQTL_mapping/batchscripts/FullGetCisSNPsJobList.sh"
    params:
        plink = "--geno --maf 0.05 --allow-no-sex --silent"
    run:
        GenesToTest = set()
        with open(input.subset_to_test, 'rU') as f:
            for line in f:
                n, gene = line.strip('\n').split('\t')
                GenesToTest.add(gene)
        fhOut = open(output.batchshellscript, 'w')
        with open(input.genes_bed, 'rU') as f:
            for i,line in enumerate (f):
                if i>=0:
                    chromosome, start, stop, gene, score, strand, *extra = line.strip('\n').split('\t')
                    if gene in GenesToTest:
                        fhOut.write("plink --bfile {} --chr {} --from-bp {} --to-bp {} --allow-extra-chr --write-snplist --out eQTL_mapping/cis_gene_snps/{} {}\n".format(input.plink_bed[:-4], chromosome, start, stop, gene, params.plink))
            fhOut.close()

rule split_GetCisSNPsJobList_toBatches:
    input:
        batchshellscript = "eQTL_mapping/batchscripts/FullGetCisSNPsJobList.sh"
    output:
        batch_out = dynamic("eQTL_mapping/batchscripts/GetCisSNPs/batch.sh.{batchid}")
    params:
        batch_size = 20
    shell:
        """
        mkdir -p eQTL_mapping/batchscripts/GetCisSNPs/
        split -{params.batch_size} {input.batchshellscript} eQTL_mapping/batchscripts/GetCisSNPs/batch.sh.
        """

rule run_GetCisSNPs_batches:
    """
    the output files of importance (the snp lists that are created by the batch
    shell script) are not included in the DAG. However, rule output points to a
    log that will only be made if batch shell script executes succesfully
    """
    input:
        batchscript = "eQTL_mapping/batchscripts/GetCisSNPs/batch.sh.{batchid}"
    output:
        "logs/eQTL_mapping/batchscripts/GetCisSNPs/batch.sh.{batchid}.log"
    shell:
        """
        mkdir -p eQTL_mapping/cis_gene_snps/
        bash {input.batchscript} && echo "done" > {output}
        """

rule collect_GetCisSNPs_batchlogs:
    input:
        dynamic("logs/eQTL_mapping/batchscripts/GetCisSNPs/batch.sh.{batchid}.log")
    output:
        "logs/eQTL_mapping/collect_GetCisSNPs_batchlogs.log"
    log:
        "logs/eQTL_mapping/collect_GetCisSNPs_batchlogs.stderr.log"
    shell:
        """
        echo 'done' > {output} 2> {log}
        """

rule GetGemmaJobList:
    input:
        subset_to_test = "eQTL_mapping/GeneListToTest.txt",
        plink_fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
        cis_snps_collected = "logs/eQTL_mapping/collect_GetCisSNPs_batchlogs.log",
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        GRM = "eQTL_mapping/Kinship/GRM.cXX.txt",
    output:
        batchshellscript = "eQTL_mapping/batchscripts/FullGemmaJobList.sh"
    params:
        gemma = ""
    run:
        with open(input.subset_to_test, 'rU') as f:
            with open(output.batchshellscript, 'w') as fout:
                for line in f:
                    n, gene = line.strip('\n').split('\t')
                    fout.write("gemma -lmm 2 -b {} -n {} -k {} -snps {} -outdir eQTL_mapping/gemma_nominal_results/ -o {} {}\n".format(input.plink_bed[:-4], n, input.GRM, "eQTL_mapping/cis_gene_snps/" + gene + ".snplist", gene, params.gemma))

rule split_GemmaJobList_toBatches:
    input:
        batchshellscript = "eQTL_mapping/batchscripts/FullGemmaJobList.sh"
    output:
        batch_out = dynamic("eQTL_mapping/batchscripts/Gemma/batch.sh.{gemma_batchid}")
    params:
        batch_size = 20
    shell:
        """
        mkdir -p eQTL_mapping/batchscripts/Gemma
        split -{params.batch_size} {input.batchshellscript} eQTL_mapping/batchscripts/Gemma/batch.sh.
        """

rule run_Gemma_batches:
    """
    the output files of importance (the snp lists that are created by the batch
    shell script) are not included in the DAG. However, rule output points to a
    log that will only be made if batch shell script executes succesfully
    """
    input:
        batchscript = "eQTL_mapping/batchscripts/Gemma/batch.sh.{gemma_batchid}"
    output:
        "logs/eQTL_mapping/batchscripts/Gemma/batch.sh.{gemma_batchid}.log"
    shell:
        """
        mkdir -p eQTL_mapping/gemma_nominal_results/
        bash {input.batchscript} && echo "done" > {output}
        """

rule collect_Gemma_batchlogs:
    input:
        dynamic("logs/eQTL_mapping/batchscripts/Gemma/batch.sh.{gemma_batchid}.log")
    output:
        "logs/eQTL_mapping/collect_Gemma_batchlogs.log"
    log:
        "logs/eQTL_mapping/collect_Gemma_batchlogs.stderr.log"
    shell:
        """
        echo 'done' > {output} 2> {log}
        """

rule MakeMatrixEQTL_pheno_table:
    input:
        plink_fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
        subset_to_test = "eQTL_mapping/GeneListToTest.txt",
    output:
        phenotypes = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt",
    log:
        "logs/eQTL_mapping/MakeMatrixEQTL_input.log"
    shell:
        """
        Rscript scripts/SubsetPhenotypesFromFamFile.R {input.plink_fam} {input.subset_to_test} {output.phenotypes} &>> {log}
        """

rule MakeMatrixEQTL_snp_table:
    input:
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
    output:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
    params:
        sed_search_delete = config["read_mapping"]["ReadGroupID_prefix"][:-1]
    shell:
        """
        plink --bfile  eQTL_mapping/plink/ForAssociationTesting --allow-extra-chr --recode A-transpose --tab --geno --memory 40000 --out eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps
        perl -lne '/^.+?\\t(.+?\\t).+?\\t.+?\\t.+?\\t.+?\\t(.+$)/ and print "$1$2" ' eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps.traw | sed '1s/_{params}//g' > {output}
        """

rule MakeMatrixEQTL_gene_loc:
    input:
        bed = "../../output/Genes.bed",
        subset_to_test = "eQTL_mapping/GeneListToTest.txt",
    output:
        gene_loc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt"
    shell:
        """
        cat {input.subset_to_test} | awk '{{ print $2 }}' | grep -w -F -f - {input.bed} | awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "gene", "chr", "start", "stop" }} {{ print $4,$1,$2,$3 }}' > {output}
        """

rule MakeMatrixEQTL_snp_loc:
    input:
        snps = "eQTL_mapping/plink/ForAssociationTesting.bim",
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
    output:
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "snp","chrom","pos" }} {{ print $2, $1,$4 }}' {input.snps} > {output.snp_locs}
        """

rule make_covariate_file:
    """
    Use R to make covariate file that matches order of fam file
    """
    input:
        EmptyFam = "../../output/ForAssociationTesting.temp.fam",
        MetadataExcel = "../../data/Metadata.xlsx",
        ExpressionMatrix = ExpressionMatrix,
        GenotypePCs = "../../output/PopulationStructure/pca.eigenvec"
    output:
        "../../output/Covariates/{NumGenotypePCs}GenotypePCs_and_{NumRNASeqPCs}RNASeqPCs.covariates"
    log:
        "logs/eQTL_mapping/make_covariate_file/{NumGenotypePCs}GenotypePCs_and_{NumRNASeqPCs}RNASeqPCs.covariates.txt.log"
    shell:
        """
        Rscript ../../analysis/20190427_MakeCovariateFiles.NoSex.R {input.EmptyFam} {input.MetadataExcel} {input.ExpressionMatrix} {wildcards.NumRNASeqPCs} {input.GenotypePCs} {wildcards.NumGenotypePCs} {output} &> {log}
        """

if config["eQTL_mapping"]["model_type"] == "lm":
    CovarianceMatrix = "eQTL_mapping/Kinship/IdentityMatrix.txt"
elif config["eQTL_mapping"]["model_type"] == "lmm":
    CovarianceMatrix = "eQTL_mapping/Kinship/GRM.cXX.txt"

rule MatrixEQTL:
    """Matrix EQTL script performs one cis-eqtl scan with real data and one
    scan with permutated sample labels for phenotypes for an empirical null."""
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        phenotypes = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt",
        gene_loc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt",
        covariates = get_covariates,
        GRM = CovarianceMatrix,
    output:
        results = "eQTL_mapping/MatrixEQTL/Results/Results.{covariate_set}.txt",
        fig = "eQTL_mapping/MatrixEQTL/Results/images/Results.{covariate_set}.png",
        permutated_fig = "eQTL_mapping/MatrixEQTL/Results/images/PermutatedResults.{covariate_set}.png",
        permuted_results = "eQTL_mapping/MatrixEQTL/Results/PermutatedResults.{covariate_set}.txt",
    log:
        "logs/eQTL_mapping/MatrixEQTL/{covariate_set}.log"
    shell:
        """
        Rscript scripts/MatrixEqtl_Cis.R {input.snps} {input.snp_locs} {input.phenotypes} {input.gene_loc} {input.covariates} {input.GRM} {output.results} {output.fig} {output.permuted_results} {output.permutated_fig} 250000 &> {log}
        """

rule MatrixEQTL_BestModelFromConfigFullResults:
    """
    Matrix eQTL with full output for every snp-gene pair. Also with best p-value for each gene.
    """
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        phenotypes = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt",
        gene_loc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt",
        covariates = config["eQTL_mapping"]["CovariatesForFullOutput"],
        GRM = CovarianceMatrix,
    output:
        results = "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt",
        BestGenePvals = "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/BestPvalueNonPermuted.txt",
        fig = "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/images/Results.png",
    log:
        "logs/eQTL_mapping/MatrixEQTL/ConfigCovariateModel.log"
    shell:
        """
        Rscript scripts/MatrixEqtl_Cis.AllPvals.R {input.snps} {input.snp_locs} {input.phenotypes} {input.gene_loc} {input.covariates} {input.GRM} {output.results} {output.fig} 250000 &> {log}
        """

rule PlotPCsVsEQTLs:
    input:
        MatrixEQTLModels
    output:
        CattedResult = "eQTL_mapping/MatrixEQTL/Results/ConcatenatedResult.txt",
        Plot = "eQTL_mapping/MatrixEQTL/Results/images/PCsVsEQTLs.pdf"
    log:
        "logs/eQTL_mapping/PlotPCsVsEQTLs.log"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' 'FNR>1 && $6<0.3 {{ print $1,$2,$3,$4,$5,$6,FILENAME  }}' {input} > {output.CattedResult}
        Rscript scripts/Plot_EQTLs_vs_PCs.R {output.CattedResult} {output.Plot}
        """

rule GetBestModel:
    input:
        CattedResult = "eQTL_mapping/MatrixEQTL/Results/ConcatenatedResult.txt",
    output:
        "eQTL_mapping/MatrixEQTL/Results/BestModelResults.txt"
    shell:
        """
        cat $(cat {input.CattedResult} | awk -F'\\t' '$6<0.1 {{print $7}}' | sort | uniq -c | sort -nr | head -1 | awk '{{print $2}}') > {output}
        """

rule GetEqtlGenotypes:
    """For checking genotype vs phenotype R-ggplot boxplots for eqtls"""
    input:
        eqtls = "eQTL_mapping/MatrixEQTL/Results/BestModelResults.txt",
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
    output:
        sig_genotypes = config["gitinclude_output"] + "MatrixEQTL_sig_genotypes.raw",
        sig_eqtls = config["gitinclude_output"] + "MatrixEQTL_sig_eqtls.txt"
    log:
        "eQTL_mapping/GetEqtlGenotypes.log"
    shell:
        """
        awk -F'\\t' 'NR==1 {{ print }} NR>1 && $6<0.1 {{ print }}' {input.eqtls} > {output.sig_eqtls}
        awk -F'\\t' 'NR>1 && $6<0.1 {{ print $1 }}' {input.eqtls} | sort | uniq -u | plink --bfile eQTL_mapping/plink/ForAssociationTesting --extract /dev/stdin --recode A --allow-extra-chr --out {config[gitinclude_output]}MatrixEQTL_sig_genotypes &> {log}
        """

rule MakeFastQTL_input:
    input:
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        gene_loc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt",
        phenotypes = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt",
    output:
        vcf = "eQTL_mapping/FastQTL/ForAssociationTesting.vcf.gz",
        vcftbi = "eQTL_mapping/FastQTL/ForAssociationTesting.vcf.gz.tbi",
        bed = "eQTL_mapping/FastQTL/ForAssociationTesting.bed.gz",
        bedtbi = "eQTL_mapping/FastQTL/ForAssociationTesting.bed.gz.tbi"
    shell:
        """
        # Make genotypes vcf
        plink2 --bfile eQTL_mapping/plink/ForAssociationTesting  --recode vcf-iid --allow-extra-chr --out eQTL_mapping/FastQTL/ForAssociationTesting
        sed -i '1 s/3$/2/' eQTL_mapping/FastQTL/ForAssociationTesting.vcf
        bgzip eQTL_mapping/FastQTL/ForAssociationTesting.vcf && tabix -p vcf {output.vcf}

        # Make phenotypes bed
        paste <(awk 'NR>1' {input.gene_loc} | sort) <(awk 'NR>1' {input.phenotypes} | sort) | perl -lne '/^.+?\\t(.+)$/ and print  "$1"' | bedtools sort -i - > eQTL_mapping/FastQTL/ForAssociationTesting.bed
        bgzip eQTL_mapping/FastQTL/ForAssociationTesting.bed && tabix -p bed {output.bed}
        """

rule eQTL_permutations:
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        phenotypes = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt",
        gene_loc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt",
        covariates = config["eQTL_mapping"]["CovariatesForFullOutput"],
        GRM = CovarianceMatrix,
    output:
        results = "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Permutations/Chunk.{n}.txt",
    params:
        InitialSeed = GetInitialSeedNumberForPermutationChunk,
        NumberPermutations = config["eQTL_mapping"]["PermutationChunkSize"]
    log:
        "logs/eQTL_mapping/MatrixEQTL/Permutations/Chunk.{n}.log"
    shell:
        """
        Rscript scripts/MatrixEqtl_Cis_Permutations.R {input.snps} {input.snp_locs} {input.phenotypes} {input.gene_loc} {input.covariates} {input.GRM} {output.results} {params.NumberPermutations} {params.InitialSeed} 250000 > {log}
        """

def GetInitialSeedNumberForPermutationChunk(wildcards):
    return int(wildcards.n) * int(config["eQTL_mapping"]["PermutationChunkSize"])

rule MergePermutationChunks:
    input:
        Chunks = expand("eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Permutations/Chunk.{n}.txt", n=range(0, int(config["eQTL_mapping"]["NumberPermutationChunks"])))
    output:
        "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/PermutationsCombined.txt"
    shell:
        """
        awk 'NR==1 && FNR=NR {{print}} FNR>1 {{input.Chunks}}' {input.Chunks} > {output}
        """

# rule prepare_vcf_for_VerifyBamID:
#     """Grab vcf from plink file for association testing, since these genotypes
#     are already well filtered for errors"""
#     input:
#         bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
#     output:
#         "qc/verifybamid/VerifyBamID.chr21.vcf.gz"
#     shell:
#         """
#         # plink2 --bfile eQTL_mapping/plink/ForAssociationTesting --chr 21 --recode vcf-fid --allow-extra-chr --out qc/verifybamid/VerifyBamID.chr21
#         # bgzip qc/verifybamid/VerifyBamID.chr21.vcf && tabix -p vcf {output}
#         """

# rule VerifyBamID:
#     """Find best sample match (among those included for eqtl calling) for each
#     RNA-seq sample (among all samples in RNAseqsamples.tsv list).
#     Note this means that if some samples are discarded from eqtl calling, this rule will pick the best"""
#     input:
#         bam = "RNASeq/STAR/{sample}/Aligned.sortedByCoord.out.bam",
#         bambai = "RNASeq/STAR/{sample}/Aligned.sortedByCoord.out.bam",
#         vcf = "filtered/21.vcf.gz"
#         # bam = "/project/gilad/bjf79/Chimp_eQTL/dna-seq-gatk-variant-calling/dedup_merged/{sample}.sorted.bam",
#         # bai = "/project/gilad/bjf79/Chimp_eQTL/dna-seq-gatk-variant-calling/dedup_merged/{sample}.sorted.bam.bai",
#         # vcf = "qc/verifybamid/VerifyBamID.chr21.vcf.gz"
#     output:
#         bestSM = "qc/verifybamid/{sample}.bestSM",
#         chr21bam = "qc/verifybamid/bams/{sample}.bam"
#     log:
#         "logs/Admixture/qc/VerifyBamID.{sample}.log"
#     shell:
#         """
#         samtools view -bh {input.bam} 21 1-1000000 > {output.chr21bam} 2> {log}
#         samtools index {output.chr21bam} &>> {log}
#         verifyBamID --vcf {input.vcf} --bam {output.chr21bam} --best --out qc/verifybamid/{wildcards.sample} &>> {log}
#         """




# TODO:
# rule permutation testing, Needs new MatrixEQTL script.
# pseudocode:
# choose number of permutations and batch size (N). Execute R script, each will output best Pvals for N permutations. Next rule to cat all those permutations. Finally, R script to get permutated gene-wise P-val.
# rule LD_decay_plot
