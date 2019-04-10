# TODO: Make covariate file (R or by hand). Make kinship matrix via gemma.
# compare to KING.Use plink2 to make bimbam file. Make phenotype matrix file
# (from count table). Use gemma for association testing ...
# Use similar strategy as leafcutter prepare juncfiles; make quick python
# wrapper script, taking as input a txt file with parameters for each gemma
# call. The parameters to include are a start and stop coordinate of gene, the
# gene name, the column number (n) in the phenotype file, the bfile, the output
# prefix. The python script will then make calls to gemma roughly like so:
#
# gemma -lmm 2 -b {bfile} -n {column} -k output/results.cXX.txt -snps <(grep {genename} eQTL_mapping/Misc/cis_gene_windows | bedtools intersect -a - -b eQTL_mapping/plink/ForAssociationTesting.snps.bed -sorted) -outdir {outdir} -o {genaname} -c {covariates}
#
# Make expression matrix for all genes to be test
# Make all bfiles/fam files for each gene to be tested.
# Divide into sets of 1000 genes to test at a time, using dynamic keyword.
# python script to make gemma parameter file

# rule all:
#     input:
#         # "eQTL_mapping/MatrixEQTL/Results.txt",
#         expand("eQTL_mapping/MatrixEQTL/Results.{PCsInModel}.txt", PCsInModel=range(0,15))
#         # "eQTL_mapping/batchscripts/FullGemmaJobList.sh",
#         # "eQTL_mapping/batchscripts/FullGetCisSNPsJobList.sh",
#         # "eQTL_mapping/GeneListToTest.txt",
#         # "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps.raw",

rule make_cis_gene_windows:
    input:
        genes_bed = config["eQTL_mapping"]["genes"]
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


rule make_plink_file_for_gemma:
    input:
        vcf = ancient("PopulationSubstructure/ReferencePanelMerged.annotated.vcf.gz"),
    output:
        bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        fam = config['gitinclude_output'] + "ForAssociationTesting.temp.fam",
        bim = "eQTL_mapping/plink/ForAssociationTesting.bim"
    log:
        "logs/eQTL_mapping/make_plink_file_for_gemma.log"
    params:
        "--keep-fam <(echo {fam})".format(fam=config["read_mapping"]["ReadGroupID_prefix"].replace('-',''))
    shell:
        """
        plink --id-delim '-' --vcf {input.vcf} --vcf-half-call m --allow-extra-chr --make-bed --out eQTL_mapping/plink/ForAssociationTesting.temp {params} &> {log}
        mv eQTL_mapping/plink/ForAssociationTesting.temp.bed {output.bed}
        mv eQTL_mapping/plink/ForAssociationTesting.temp.bim {output.bim}
        mv eQTL_mapping/plink/ForAssociationTesting.temp.fam {output.fam}
        """

rule make_filtered_phenotype_matrix:
    """
    Use R to filter TPM count matrix for genes to map eqtls for.
    """
    input:
        fam = config['gitinclude_output'] + "ForAssociationTesting.temp.fam",
        counttable = config["gitinclude_output"] + "CountTable.tpm.txt.gz",
        genes_bed = config["eQTL_mapping"]["genes"]
    output:
        fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
        gene_list = "eQTL_mapping/plink/ForAssociationTesting.fam.phenotype-list"
    log:
        "logs/eQTL_mapping/make_filtered_phenotype_matrix.log"
    shell:
        """
        Rscript ../../analysis/20190327_MakeFamPhenotypeFile.R {input.counttable} {input.fam} {output.fam} {output.gene_list} {input.genes_bed} &> {log}
        """

rule Make_GRM:
    """Genetetic relatedness matrix for gemma LMM"""
    input:
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
    output:
        GRM = "eQTL_mapping/GRM.cXX.txt"
    log:
        "logs/eQTL_mapping/make_GRM.log"
    shell:
        """
        gemma -gk 1 -bfile eQTL_mapping/plink/ForAssociationTesting -o GRM -outdir eQTL_mapping/ &> {log}
        """

rule Make_King_GRM:
    """Genetic relatedness matrix from plink2 (king robust)"""
    input:
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
    output:
        GRM = "eQTL_mapping/GRM.king.txt"
    log:
        "logs/eQTL_mapping/make_King_GRM.log"
    shell:
        """
        plink2 --bfile eQTL_mapping/plink/ForAssociationTesting --king
        """

# rule make_covariate_file:
#     """
#     Use R to make covariate file that matches order of fam file
#     """
#     input:
#     output:
#     log:
#     shell:

#put subset gene list here. to subset fam.phenotype-list instead of gemma command batch script. That way MatrixEQTL isn't dependent on the GetCisSnp steps.

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
        GRM = "eQTL_mapping/GRM.cXX.txt",
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
        plink --maf 0.05 --bfile  eQTL_mapping/plink/ForAssociationTesting --allow-extra-chr --recode A-transpose --tab --geno --memory 40000 --out eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps
        perl -lne '/^.+?\\t(.+?\\t).+?\\t.+?\\t.+?\\t.+?\\t(.+$)/ and print "$1$2" ' eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps.traw | sed '1s/_{params}//g' > {output}
        """

rule MakeMatrixEQTL_gene_loc:
    input:
        bed = config["eQTL_mapping"]["genes"],
        subset_to_test = "eQTL_mapping/GeneListToTest.txt",
    output:
        gene_loc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt"
    shell:
        """
        cat {input.subset_to_test} | awk '{{ print $2 }}' | grep -w -F -f - {input.bed} | awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "gene", "chr", "start", "stop" }} {{ print $4,$1,$2,$3 }}' > {output}
        """

rule MakeMatrixEQTL_snp_loc:
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
    output:
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' 'NR==1 {{ print "snp", "chrom", "pos" }} NR>1 {{ split($1,a,":"); print $1, a[1], a[2] }}' {input.snps} > {output.snp_locs}
        """

rule MatrixEQTL:
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        phenotypes = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt",
        gene_loc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt",
        covariates = "../../output/Covariates.{PCsInModel}.txt",
        grm = "scratch/plink2.king.Reformatted",
    output:
        results = "eQTL_mapping/MatrixEQTL/Results/QN.{PCsInModel}.txt",
        fig = "eQTL_mapping/MatrixEQTL/Results/QN.{PCsInModel}.png"
    log:
        "logs/eQTL_mapping/MatrixEQTL.{PCsInModel}.log"
    shell:
        """
        Rscript scripts/MatrixEqtl_Cis.R {input.snps} {input.snp_locs} {input.phenotypes} {input.gene_loc} {input.covariates} {input.grm} {output.results} {output.fig} &> {log}
        """

# rule which result maximizes eqtls or egenes
# rule permutation testing
# rule 
