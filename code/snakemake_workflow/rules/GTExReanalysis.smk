# fastQTL --vcf genotypes.vcf.gz --bed phenotypes.bed.gz --region 22:17000000-18000000 --permute 1000 --out permutations.default.txt.gzfastQTL --vcf genotypes.vcf.gz --bed phenotypes.bed.gz --region 22:17000000-18000000 --permute 1000 --out permutations.default.txt.gz
rule MakeSubsamples_VaryingSize:
    input: config["GTEx_eQTL_mapping"]["sample_list"]
    output:
        "GTEX_renalysis/SubsamplesVaryingSize/SampleLists/{n}.txt"
    shell:
        """
        shuf -n {wildcards.n} {input} > {output}
        """

rule MakeVCF_forSubsample:
    """
    FastQTL throws an error when there is no variation for genotype or phenotypes.
    When doing subsample, it is possible that some genotypes will have no
    variation. The solution implemented in this rule is to create a new vcf for
    each subsample with a MAF threshold based on the subsample.
    """
    input:
        vcf = config["GTEx_eQTL_mapping"]["vcf"],
        vcf_tbi = config["GTEx_eQTL_mapping"]["vcf"] + ".tbi",
        samples = "GTEX_renalysis/SubsamplesVaryingSize/SampleLists/{n}.txt"
    output:
        vcf = "GTEX_renalysis/SubsamplesVaryingSize/vcf/{n}.vcf.gz",
        tbi =  "GTEX_renalysis/SubsamplesVaryingSize/vcf/{n}.vcf.gz.tbi"
    log:
        "logs/GTEX_renalysis/MakeVCF_forSubsample/{n}.log"
    shell:
        """
        bcftools view -O z -S {input.samples} -q 0.1:minor {input.vcf} > {output.vcf} 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

def GetNumberPEER_factors(wildcards):
    N = int(wildcards.n)
    if N<50:
        return(10)
    elif N<150:
        return(15)
    elif N<250:
        return(30)

rule SubsetCovariatesFiles:
    """
    Similar to recommendations of GTEx, use less covariates for smaller sample
    sizes. Save subset of covariates to new file.
    """
    input:
        covariates = config["GTEx_eQTL_mapping"]["covariates"]
    params:
        GetNumberPEER_factors
    output:
        "GTEX_renalysis/SubsamplesVaryingSize/covariates/{n}.txt"
    log: "logs/GTEX_renalysis/SubsetCovariatesFiles/{n}.log"
    shell:
        """
        set +o pipefail;
        grep -v "InferredCov" {input.covariates} > {output} 2> {log}
        grep "InferredCov" {input.covariates} | head -n {params} >> {output} 2>> {log}
        """

rule FastQTL_GTEx_varyingSampleSize:
    input:
        samples = "GTEX_renalysis/SubsamplesVaryingSize/SampleLists/{n}.txt",
        vcf = "GTEX_renalysis/SubsamplesVaryingSize/vcf/{n}.vcf.gz",
        vcf_tbi = "GTEX_renalysis/SubsamplesVaryingSize/vcf/{n}.vcf.gz.tbi",
        phenotypes = config["GTEx_eQTL_mapping"]["phenotypes"],
        phenotypes_tbi = config["GTEx_eQTL_mapping"]["phenotypes"] + ".tbi",
        covariates = "GTEX_renalysis/SubsamplesVaryingSize/covariates/{n}.txt"
    output:
        permutation_chunk = "GTEX_renalysis/SubsamplesVaryingSize/output_chunks/SampleSize_{n}.Chunk_{chunk}.txt.gz"
    log:
        "logs/GTEx_reanalysis/FastQTL_GTEx/{n}.{chunk}.log"
    params:
        config["GTEx_eQTL_mapping"]["FastQTL_params"]
    shell:
        """
        ~/software/FastQTL/bin/fastQTL.static --cov {input.covariates} --vcf {input.vcf} --include-samples {input.samples} --bed {input.phenotypes} --chunk {wildcards.chunk} {config[GTEx_eQTL_mapping][chunks]} --permute 1000 10000 {params} --out {output.permutation_chunk} &> {log}
        """

rule MergeFastQTL_GTEx_varyingSampleSize:
    input:
        permutation_chunk = expand("GTEX_renalysis/SubsamplesVaryingSize/output_chunks/SampleSize_{{n}}.Chunk_{chunk}.txt.gz", chunk=range(1, 1+int(config["GTEx_eQTL_mapping"]["chunks"])))
    output:
        "GTEX_renalysis/SubsamplesVaryingSize/output/SampleSize_{n}.txt.gz"
    shell:
        """
        cat {input} > {output}
        """

rule AddQvalueColumn:
    input:
        "GTEX_renalysis/SubsamplesVaryingSize/output/SampleSize_{n}.txt.gz"
    output:
        "../../output/GTEX_renalysis/SampleSize_{n}.txt.gz"
    shell:
        """
        Rscript scripts/AddQvaluesToFastQTLPermutationOutput.R {input} ../../output/GTEX_renalysis/SampleSize_{wildcards.n}.txt
        gzip ../../output/GTEX_renalysis/SampleSize_{wildcards.n}.txt
        """

