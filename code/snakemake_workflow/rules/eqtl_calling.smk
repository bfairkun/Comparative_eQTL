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

rule all:
    input:
        "eQTL_mapping/batchscripts/FullGemmaJobList.sh"

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
    output:
        fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
        gene_list = "eQTL_mapping/plink/ForAssociationTesting.fam.phenotype-list"
    log:
        "logs/eQTL_mapping/make_filtered_phenotype_matrix.log"
    shell:
        """
        Rscript ../../analysis/20190327_MakeFamPhenotypeFile.R {input.counttable} {input.fam} {output.fam} {output.gene_list} &> {log}
        """

# rule make_covariate_file:
#     """
#     Use R to make covariate file that matches order of fam file
#     """
#     input:
#     output:
#     log:
#     shell:

rule GetCisSNPsJobList:
    input:
        genes_bed = "eQTL_mapping/Misc/cis_gene_windows.bed",
        genes_to_test = "eQTL_mapping/plink/ForAssociationTesting.fam.phenotype-list",
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
    output:
        batchshellscript = "eQTL_mapping/batchscripts/FullGetCisSNPsJobList.sh"
    params:
        plink = "--geno 0.05 --maf 0.05 --allow-no-sex --silent"
    run:
        with open(input.genes_to_test, 'rU') as f:
            GenesToTest = set(f.read().splitlines())
        fhOut = open(output.batchshellscript, 'w')
        with open(input.genes_bed, 'rU') as f:
            for i,line in enumerate (f):
                if i:
                    chromosome, start, stop, gene, score, strand, *extra = line.strip('\n').split('\t')
                    if gene in GenesToTest:
                        fhOut.write("plink --bfile {} --chr {} --from-bp {} --to-bp {} --allow-extra-chr --write-snplist --out eQTL_mapping/cis_gene_snps/{} {}\n".format(input.plink_bed[:-4], chromosome, start, stop, gene, params.plink))
            fhOut.close()

rule split_GetCisSNPsJobList_toBatches:
    input:
        batchshellscript = "eQTL_mapping/batchscripts/FullGetCisSNPsJobList.sh"
    output:
        batch_out = dynamic("eQTL_mapping/batchscripts/GetCisSNPs.batch.sh.{batchid}")
    params:
        batch_size = 200
    shell:
        """
        split -{params.batch_size} {input.batchshellscript} eQTL_mapping/batchscripts/GetCisSNPs.batch.sh.
        """

rule run_GetCisSNPs_batches:
    """
    the output files of importance (the snp lists that are created by the batch
    shell script) are not included in the DAG. However, rule output points to a
    log that will only be made if batch shell script executes succesfully
    """
    input:
        batchscript = "eQTL_mapping/batchscripts/GetCisSNPs.batch.sh.{batchid}"
    output:
        "logs/eQTL_mapping/batchscripts/GetCisSNPs.batch.sh.{batchid}.log"
    shell:
        """
        mkdir -p eQTL_mapping/cis_gene_snps/
        bash {input.batchscript} && echo "done" > {output}
        """

rule collect_GetCisSNPs_batchlogs:
    input:
        dynamic("logs/eQTL_mapping/batchscripts/GetCisSNPs.batch.sh.{batchid}.log")
    output:
        "logs/eQTL_mapping/collect_GetCisSNPs_batchlogs.log"
    shell:
        """
        echo 'done' > {output}
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
        gemma -gk 1 -bfile {plink_bed} -o GRM -outdir eQTL_mapping/ &> {log}
        """

rule GetGemmaJobList:
    input:
        genes_to_test = "eQTL_mapping/plink/ForAssociationTesting.fam.phenotype-list",
        plink_fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
        # cis_snps_collected = "logs/eQTL_mapping/collect_GetCisSNPs_batchlogs.log"
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        GRM = "eQTL_mapping/GRM.cXX.txt",
    output:
        batchshellscript = "eQTL_mapping/batchscripts/FullGemmaJobList.sh"
    params:
        gemma = ""
    run:
        with open(input.genes_to_test, 'rU') as f:
            with open(output.batchshellscript, 'w') as fout:
                for i,line in enumerate (f):
                    n = i+1
                    gene = line.strip('\n')
                    if i < 10:
                        fout.write("gemma -lmm 2 -b {} -n {} -k {} -snps {} -outdir eQTL_mapping/gemma_nominal_results/ -o {} {}\n".format(input.plink_bed[:-4], n, input.GRM, "eQTL_mapping/cis_gene_snps/" + gene + ".snplist", gene, params.gemma))


# gemma -lmm 2 -b {bfile} -n {column} -k output/results.cXX.txt -snps <(grep {genename} eQTL_mapping/Misc/cis_gene_windows | bedtools intersect -a - -b eQTL_mapping/plink/ForAssociationTesting.snps.bed -sorted) -outdir {outdir} -o {genaname} -c {covariates}

