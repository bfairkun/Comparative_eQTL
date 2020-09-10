rule SubsetVcfForCellCompositionGWAS:
    input:
        vcf = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz",
        tbi = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz.tbi",
        pheno = "../../output/CellProportionPhenotypesNormalizedForGWAS.tab"
    output:
        vcf = "AddressReviews/GWAS/GTEX.v8.vcf.gz"
    shell:
        """
        bcftools view --force-samples -S <(awk -F'\\t' 'NR>1 {{print $1}}' {input.pheno} | sort | uniq) -O z -m 2 -M 2 -q 0.05:minor {input.vcf} > {output.vcf}
        """

rule MakePlinkBedForGWAS:
    input:
        vcf = "AddressReviews/GWAS/GTEX.v8.vcf.gz"
    output:
        bed =  "AddressReviews/plink/GTEX.v8.bed",
        fam = "AddressReviews/plink/GTEX.v8.fam",
        bim = "AddressReviews/plink/GTEX.v8.bim"
    log:
        "logs/plink/GWASMakeBed"
    shell:
        """
        plink --vcf {input.vcf} --vcf-half-call m  --hwe 1e-7.5 --make-bed --id-delim "-" --out AddressReviews/plink/GTEX.v8 &> {log}
        """

rule MakeCovariatesForGWAS:
    input:
        bed =  "AddressReviews/plink/GTEX.v8.bed",
        fam = "AddressReviews/plink/GTEX.v8.fam",
        bim = "AddressReviews/plink/GTEX.v8.bim",
        covariates = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Heart_Left_Ventricle.v8.covariates.txt",
        pheno = "../../output/CellProportionPhenotypesNormalizedForGWAS.tab"
    output:
        bed =  "AddressReviews/GWAS/GTEX.v8.bed",
        fam = "AddressReviews/GWAS/GTEX.v8.fam",
        bim = "AddressReviews/GWAS/GTEX.v8.bim",
        covariates = "AddressReviews/GWAS/GTEX.v8.cov",
    shell:
        """
        cp {input.bed} {output.bed}
        cp {input.bim} {output.bim}
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/CovariatesForCellCompositionGWAS.R
        """

rule GemmaForCellTypeQTL:
    input:
        bed =  "AddressReviews/GWAS/GTEX.v8.bed",
        fam = "AddressReviews/GWAS/GTEX.v8.fam",
        bim = "AddressReviews/GWAS/GTEX.v8.bim",
        covariates = "AddressReviews/GWAS/GTEX.v8.cov",
    output:
        "gemma out"
    shell:
        """
        gemma -lm 1 -bfile AddressReviews/GWAS/GTEX.v8 -c {input.covariates}
        """
