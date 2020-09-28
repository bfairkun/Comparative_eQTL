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
        covariates_lmm = "AddressReviews/GWAS/GTEX.v8.lmm.cov",
    shell:
        """
        cp {input.bed} {output.bed}
        cp {input.bim} {output.bim}
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/CovariatesForCellCompositionGWAS.R
        awk -F'\\t' '{{print $1, $NF}}' {output.covariates} > {output.covariates_lmm}
        """

rule GemmaForCellTypeGRM:
    input:
        bed =  "AddressReviews/GWAS/GTEX.v8.bed",
        fam = "AddressReviews/GWAS/GTEX.v8.fam",
        bim = "AddressReviews/GWAS/GTEX.v8.bim",
    output:
        "output/result.cXX.txt"
    shell:
        "gemma -bfile AddressReviews/GWAS/GTEX.v8 -gk 1"

rule GemmaForCellTypeQTL:
    input:
        bed =  "AddressReviews/GWAS/GTEX.v8.bed",
        fam = "AddressReviews/GWAS/GTEX.v8.fam",
        bim = "AddressReviews/GWAS/GTEX.v8.bim",
        covariates = "AddressReviews/GWAS/GTEX.v8.lmm.cov",
        grm = "output/result.cXX.txt"
    output:
        "output/result.assoc.txt.gz"
    shell:
        """
        gemma -lmm 1 -bfile AddressReviews/GWAS/GTEX.v8 -c {input.covariates} -k {input.grm} -km 1 
        gzip output/result.assoc.txt
        """

rule GetGTExGenesAndDispersionGenes_HumanBed:
    input:
        GtexGenes = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Heart_Left_Ventricle.v8.normalized_expression.bed.gz",
        DispersionGenes = "../../output/OverdispersionEstimatesFromChimp.NoLengthNorm.txt"
    output:
        "Misc/GtexDispersionGenes.hg38.bed"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' 'NR>1 {{ print $1}}' {input.DispersionGenes} | grep -f - <(zcat {input.GtexGenes} | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3,$4}}') | bedtools sort -i - > {output} 
        """

rule GetRawCountTableToOutput:
    input:
        "RNASeq/STAR/CountTable.txt.gz"
    output:
        "../../output/STAR.RawCountTable.txt.gz"
    shell:
        "zcat {input} | awk -v OFS='\\t' '{{print $0}}' | gzip - > {output}"

rule BedtoolsClosestGeneToGWASSNPs:
    input:
        Control = "../../output/CellProportionGWAS.RandomControlloci.bed",
        Test = "../../output/CellProportionGWAS.loci.bed",
        bed = "Misc/GtexDispersionGenes.hg38.bed"
    output:
        Control = "../../output/CellProportionGWAS.RandomControlloci.closestGenes.txt",
        Test = "../../output/CellProportionGWAS.Testloci.closestGenes.txt"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '{{print "chr"$1, $2,$3,$4,$5,$6}}' {input.Control} | bedtools sort -i - | bedtools closest -a -  -b {input.bed} -d | awk -F'\\t' '$NF<100000 {{ print $(NF-1), $NF }}' > {output.Control}
        awk -F'\\t' -v OFS='\\t' '{{print "chr"$1, $2,$3,$4,$5,$6}}' {input.Test} | bedtools sort -i - | bedtools closest -a -  -b {input.bed} -d | awk -F'\\t' '$NF<100000 {{ print $(NF-1), $NF }}' > {output.Test}
        """

rule BootstrapErrorForTEDDispersion:
    input:
        "MiscOutput/NormalizedExpressionPerCellType.rds"
    output:
        "AddressReviews/TED_Bootstraps/{seed}.rds"
    log:
        "logs/AddressReviews/BootstrapError.{seed}.log"
    params:
        ChunkSize = DispersionTEDBootstrapChunkSize
    shell:
        """
        /software/R-3.5.1-el7-x86_64/bin/Rscript scripts/BootstrapErrorDispersionTED.R {wildcards.seed} {params.ChunkSize}
        """

rule CombineBootstrapRepsForTEDDispersion:
    input:
        expand("AddressReviews/TED_Bootstraps/{seed}.rds", seed=DispersionTEDBootstrapInitialSeeds),
    output:
        "../../output/CellTypeDispersion.SE.tsv.gz"
    shell:
        """
        /software/R-3.5.1-el7-x86_64/bin/Rscript scripts/CombineTEDBootstrapReps.R
        """
