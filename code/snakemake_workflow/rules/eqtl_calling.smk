# TODO: Make covariate file (R or by hand). Make kinship matrix via gemma.
# compare to KING.Use plink2 to make bimbam file. Make phenotype matrix file
# (from count table). Use gemma for association testing ...
# Use similar strategy as leafcutter prepare juncfiles; make quick python
# wrapper script, taking as input a txt file with parameters for each gemma
# call. The parameters to include are a start and stop coordinate of gene, the
# gene name, the column number (n) in the phenotype file, the bfile, the output
# prefix. The python script will then make calls to gemma roughly like so:
#
# gemma -lmm 2 -b {bfile} -n {column} -k output/results.cXX.txt -snps <(bedtools intersect -a <(bedtools flank -l 1000000 -r 1000000 -s) -b SnpList.bed -sorted -g genome.faidx | awk -F'\t' '{print $4}') -outdir {outdir} -o {genaname} -c {covariates}
#
# Make expression matrix for all genes to be test
# Make all bfiles/fam files for each gene to be tested.
# Divide into sets of 1000 genes to test at a time, using dynamic keyword.
# python script to make gemma parameter file

rule make_cis_gene_windows:
    input:
        config["eQTL_mapping"]["gene_bed"]
    output:
        "eQTL_mapping/Misc/cis_gene_windows.bed"
    log:
        "logs/eQTL_mapping/cis_gene_windows.log"
    params:
        MaxDistFromTSS = config["eQTL_mapping"]["TSS_flanking_snp_width"],
        faidx = config["ref"]["genome"] + ".fai",
    shell:
        """
        bedtools sort -i {input} | bedtools flank -l 1 -r 0 -s -g {params.faidx} | bedtools slop -b {params.MaxDistFromTSS} -g {params.faidx} > {output}
        """

rule make_snp_bed:
    "UCSC format bed, snp-id is 4th column"
    input:
        bim = "eQTL_mapping/plink/ForAssociationTesting.bim"
    output:
        "eQTL_mapping/plink/ForAssociationTesting.snps.bed"
    log:
        "logs/eQTL_mapping/make_snp_bed.log"
    shell:
        """
        awk -F'\\t' '{{ print $1,$4,$4+1,$2,".", "." }}'
        """


rule make_plink_file_for_gemma:
    input:
        vcf = ancient("PopulationSubstructure/ReferencePanelMerged.annotated.vcf.gz"),
    output:
        bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        fam = "eQTL_mapping/plink/ForAssociationTesting.temp.fam",
        bim = "eQTL_mapping/plink/ForAssociationTesting.bim"
    log:
        "logs/eQTL_mapping/make_plink_file_for_gemma.log"
    params:
        MyParams = "--keep-fam <(echo {fam})".format(fam=config["read_mapping"]["ReadGroupID_prefix"])
    shell:
        """
        plink --vcf {input.vcf} --vcf-half-call m --allow-extra-chr --make-bed --out eQTL_mapping/plink/ForAssociationTesting.temp {params.MyParams} &> {log}}
        mv eQTL_mapping/plink/ForAssociationTesting.temp.bed eQTL_mapping/plink/ForAssociationTesting.bed
        mv eQTL_mapping/plink/ForAssociationTesting.temp.bim eQTL_mapping/plink/ForAssociationTesting.bim
        """

rule make_filtered_phenotype_matrix:
    """
    Use R to filter TPM count matrix for genes to map eqtls for.
    """
    input:
    output:
        ".fam"
    log:
    shell:

rule make_covariate_file:
    """
    Use R to make covariate file that matches order of fam file
    """
    input:
    output:
    log:
    shell:

rule make_gemma_wrapper_parameter_file:
    input:
        cis_gene_bed = "",
    output:
        dynamic()
    log:
    shell:
        """
        python
        """

rule gemma_wrapper:
    input:
        parameter_file = ""
    output:
        "eQTL-mapping/Associations/logs/{}"
    shell:
        """
        shell
        """

