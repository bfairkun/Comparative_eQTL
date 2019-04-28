rule all:
    input:
        "../../output/Genes.bed",
        "../../output/log10TPM.StandardizedAndNormalized.txt",
        "eQTL_mapping/MatrixEQTL/Results.BestModelResults.txt"

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

rule make_plink_file_for_testing:
    """
    .fam file output to the gitinclude output directory for convenience, so
    that I can easily access it from my local laptop with git, to edit scripts
    that choose covariates and such with my local RStudio... For the covariate
    files must have samples ordered same as the fam (phenotype) file
    """
    input:
        vcf = ancient("PopulationSubstructure/ReferencePanelMerged.annotated.vcf.gz"),
    output:
        bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        fam = config['gitinclude_output'] + "ForAssociationTesting.temp.fam",
        bim = "eQTL_mapping/plink/ForAssociationTesting.bim"
    log:
        "logs/eQTL_mapping/make_plink_file_for_gemma.log"
    params:
        stdparams = "--keep-fam <(echo {fam})".format(fam=config["read_mapping"]["ReadGroupID_prefix"].replace('-','')),
        extra = "--remove <(printf Pan_troglodytes_ThisStudy\\tMD_And\\n) --maf 0.05 --geno --hwe 1e-7.5"
    shell:
        """
        plink --id-delim '-' --vcf {input.vcf} --vcf-half-call m --allow-extra-chr --make-bed --out eQTL_mapping/plink/ForAssociationTesting.temp {params.stdparams} {params.extra} &> {log}
        mv eQTL_mapping/plink/ForAssociationTesting.temp.bed {output.bed}
        mv eQTL_mapping/plink/ForAssociationTesting.temp.bim {output.bim}
        mv eQTL_mapping/plink/ForAssociationTesting.temp.fam {output.fam}
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

rule make_filtered_phenotype_matrix:
    """
    Use R to filter TPM count matrix for genes to map eqtls for. Additionally,
    do quantile normalization here.
    """
    input:
        fam = config['gitinclude_output'] + "ForAssociationTesting.temp.fam",
        counttable = config["gitinclude_output"] + "CountTable.tpm.txt.gz",
        genes_bed = "../../output/Genes.bed",
        transcripts_to_genes = "../../data/Biomart_export.Pan_Tro_3.geneids.txt"
    output:
        log10TPM = "../../output/log10TPM.txt",
        standardized = "../../output/log10TPM.StandardizedAndNormalized.txt",
        fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
        gene_list = "eQTL_mapping/plink/ForAssociationTesting.fam.phenotype-list"
    log:
        "logs/eQTL_mapping/make_filtered_phenotype_matrix.log"
    shell:
        """
        # Rscript to filter for autosomal genes and sum transcript TPM to gene TPM
        Rscript ../../analysis/20190327_MakeFamPhenotypeFile.R {input.counttable} {input.fam} {output.log10TPM} {input.genes_bed} {input.transcripts_to_genes} &> {log}

        # script to standardize and quantile normalize
        python3 scripts/StandardizeAndQuantileNormalize.py {output.log10TPM} {output.standardized} &>> {log}

        # convert to .fam format expected by gemma
        Rscript scripts/MergeFamWithPhenotypes.R {output.standardized} {input.fam} {output.fam} {output.gene_list} &>> {log}
        """

rule prune_plink_files_for_GRM:
    input:
        bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
    output:
        bed =  "eQTL_mapping/Kinship/ForAssocationTesting.pruned.bed",
    log:
        "logs/eQTL_mapping/prune_plink_files.log"
    shell:
        """
        plink --bfile eQTL_mapping/plink/ForAssociationTesting --allow-extra-chr --indep-pairwise 50 5 0.5 &> {log}
        plink --bfile eQTL_mapping/plink/ForAssociationTesting --allow-extra-chr --extract plink.prune.in --make-bed --out eQTL_mapping/Kinship/ForAssociationTesting.pruned &> {log}
        rm plink.prune.in plink.prune.out
        """

rule Make_GRM:
    """Genetetic relatedness matrix for gemma LMM"""
    input:
        bed =  "eQTL_mapping/Kinship/ForAssocationTesting.pruned.bed",
    output:
        GRM = "eQTL_mapping/Kinship/GRM.cXX.txt"
    log:
        "logs/eQTL_mapping/make_GRM.log"
    shell:
        """
        gemma -gk 1 -bfile eQTL_mapping/Kinship/ForAssociationTesting -o GRM -outdir eQTL_mapping/Kinship &> {log}
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
        ExpressionMatrix = "../../output/log10TPM.StandardizedAndNormalized.txt",
        GenotypePCs = "../../output/PopulationStructure/pca.eigenvec"
    output:
        "../../output/Covariates/{NumGenotypePCs}GenotypePCs_and_{NumRNASeqPCs}RNASeqPCs.covariates"
    log:
        "logs/eQTL_mapping/make_covariate_file/{NumGenotypePCs}GenotypePCs_and_{NumRNASeqPCs}RNASeqPCs.covariates.txt.log"
    shell:
        """
        Rscript ../../analysis/20190427_MakeCovariateFiles.R {input.EmptyFam} {input.MetadataExcel} {input.ExpressionMatrix} {wildcards.NumRNASeqPCs} {input.GenotypePCs} {wildcards.NumGenotypePCs} {output} &> {log}
        """


rule MatrixEQTL:
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        phenotypes = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt",
        gene_loc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt",
        covariates = get_covariates,
        # GRM = "eQTL_mapping/Kinship/IdentityMatrix.txt",
        GRM = "eQTL_mapping/Kinship/GRM.cXX.txt",
    output:
        results = "eQTL_mapping/MatrixEQTL/Results/Results.{covariate_set}.txt",
        fig = "eQTL_mapping/MatrixEQTL/Results/Results.{covariate_set}.png"
    log:
        "logs/eQTL_mapping/MatrixEQTL/{covariate_set}.log"
    shell:
        """
        Rscript scripts/MatrixEqtl_Cis.R {input.snps} {input.snp_locs} {input.phenotypes} {input.gene_loc} {input.covariates} {input.GRM} {output.results} {output.fig} &> {log}
        """

rule PickBestMatrixEQTLModelResults:
    input:
        MatrixEQTLModels
    output:
        FullResults = "eQTL_mapping/MatrixEQTL/Results.BestModelResults.txt"
    log:
        "logs/eQTL_mapping/MatrixEQTL/PickBestMatrixEQTLModelResults.txt"
    run:
        from shutil import copyfile
        from collections import Counter
        NumEqtls = Counter()
        print ( list(input) )
        for filepath in list(input):
            with open( filepath, 'rU' ) as f:
                for i,line in enumerate(f):
                    if i>=1:
                        snp, gene, beta, tstsat, p, fdr = line.strip('\t').split('\t')
                        if float(fdr) < 0.1:
                            print('yes')
                            NumEqtls[filepath] += 1
        # print( NumEqtls )
        with open(log[0], 'w') as f:
            for filepath, eqtl_count in NumEqtls.items():
                f.write("{}\t{}\n".format( filepath, eqtl_count ))
        filepath_w_most_eqtls = max(NumEqtls, key=lambda filepath: NumEqtls[filepath])
        copyfile( filepath_w_most_eqtls, str(output.FullResults) )

rule GetEqtlGenotypes:
    """For checking genotype vs phenotype R-ggplot boxplots for eqtls"""
    input:
        eqtls = "eQTL_mapping/MatrixEQTL/Results.BestModelResults.txt",
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
        plink2 --bfile {input.plink_bed} --recode vcf-fid --allow-extra-chr --out eQTL_mapping/FastQTL/ForAssociationTesting
        bgzip eQTL_mapping/FastQTL/ForAssociationTesting.vcf && tabix -p vcf {output.vcf}

        # Make phenotypes bed
        paste <(awk 'NR>1' {input.gene_loc} | sort) <(awk 'NR>1' {input.phenotypes} | sort) | perl -lne '/^.+?\\t(.+)$/ and print  "$1"' | bedtools sort -i - > eQTL_mapping/FastQTL/ForAssociationTesting.bed
        bgzip eQTL_mapping/ForAssociationTesting.bed && tabix -p bed {output.bed}
        """



# TODO:
# rule permutation testing, Needs new MatrixEQTL script.
# pseudocode:
# choose number of permutations and batch size (N). Execute R script, each will output best Pvals for N permutations. Next rule to cat all those permutations. Finally, R script to get permutated gene-wise P-val.
# rule LD_decay_plot
