# fastQTL --vcf genotypes.vcf.gz --bed phenotypes.bed.gz --region 22:17000000-18000000 --permute 1000 --out permutations.default.txt.gzfastQTL --vcf genotypes.vcf.gz --bed phenotypes.bed.gz --region 22:17000000-18000000 --permute 1000 --out permutations.default.txt.gz

rule DownloadGtexData:
    """
    TODO. finish this wget rule for reproducibility sake. In the meantime I
    just manually executed this rule but did not include it in this snakemake
    pipeline
    """
    output:
        expand("GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/{tissue}.v8.normalized_expression.bed.gz.tbi", tissue=GTExTissues),
        expand("GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/{tissue}.v8.covariates.txt", tissue=GTExTissues)
    shell:
        """
        #incomplete
        wget
        """

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


rule SubsetTSP_AndMatchedSnps:
    input:
        MatchedSnps = "GTEX_renalysis/LeflerSnps/ControlMatrixEQTL/ControlSnps.loc",
        LeflerSnps = "../../data/LeflerSharedPolymorphisms.hg38.bed",
        GtexVcf = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz",
        GtexVcfTbi = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz.tbi",
    output:
        GenotypeRaw = "GTEX_renalysis/LeflerSnps/Genotypes.Control.vcf.gz",
        tbi = "GTEX_renalysis/LeflerSnps/Genotypes.Control.vcf.gz.tbi",
        GenotypeRawLeffler = "GTEX_renalysis/LeflerSnps/Genotypes.TSP.vcf.gz",
        tbiLeffler = "GTEX_renalysis/LeflerSnps/Genotypes.TSP.vcf.gz.tbi",
    shell:
        """
        bcftools view -q 0.01:minor -R <(cat {input.MatchedSnps} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{ print "chr"$2, $3, $3+1 }}') {input.GtexVcf} | bcftools sort -O z > {output.GenotypeRaw}
        tabix -p vcf {output.GenotypeRaw}
        bcftools view -q 0.01:minor -R <(cat {input.LeflerSnps}) {input.GtexVcf} | bcftools sort -O z > {output.GenotypeRawLeffler}
        tabix -p vcf {output.GenotypeRawLeffler}
        """

rule FastQTL_GTEx_TSP:
    input:
        Vcf = "GTEX_renalysis/LeflerSnps/Genotypes.{SnpSet}.vcf.gz",
        tbi = "GTEX_renalysis/LeflerSnps/Genotypes.{SnpSet}.vcf.gz.tbi",
        bed = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/{tissue}.v8.normalized_expression.bed.gz",
        bedtbi = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/{tissue}.v8.normalized_expression.bed.gz.tbi",
        covariates = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/{tissue}.v8.covariates.txt"
    log:
        "logs/FastQTL_GTEx_TSP.{tissue}.{SnpSet}.log"
    output:
        "eQTL_mapping/SharedPolymorphisms/GTExAllTissues/{tissue}.{SnpSet}.txt.gz"
    wildcard_constraints:
        SnpSet = "TSP|Control",
        tissue = "|".join(GTExTissues)
    shell:
        """
        mkdir -p {config[temp_files_prefix]}
        set +e
        # fastqtl exits with error if no variants in chunk found. we want to allow that
        for i in $(seq 1 30); do
            ~/software/FastQTL/bin/fastQTL.static --cov {input.covariates} --vcf {input.Vcf}  --bed {input.bed} -L {log}  --chunk $i 30 --out {config[temp_files_prefix]}{wildcards.tissue}.{wildcards.SnpSet}.$i
        done
        set -e
        cat {config[temp_files_prefix]}{wildcards.tissue}.{wildcards.SnpSet}.*  | gzip -c > {output}
        """

rule Merge_FastQTL_TSP:
    input:
        expand("eQTL_mapping/SharedPolymorphisms/GTExAllTissues/{tissue}.{SnpSet}.txt.gz", tissue=GTExTissues, SnpSet=['TSP', 'Control']),
    output:
        "eQTL_mapping/SharedPolymorphisms/GTExAllTissues.Combined.txt"
    run:
        import gzip
        with open(output[0], "w") as out:
            for f in input:
                with gzip.open(f, 'rt') as f_in:
                    for line in f_in:
                        out.write(' '.join([f] + [line] ))

rule FastQTL_TSP_compress:
    input:
        "eQTL_mapping/SharedPolymorphisms/GTExAllTissues.Combined.txt"
    output:
        "../../output/TSP.eQTLs.GTEx.AllTissues.txt.gz"
    shell:
        """
        cat {input} | gzip - > {output}
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

rule QQPlot:
    input:
        "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt",
        "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.Permuted.txt",
    output:
        "../../output/QQPlot.png"
    log:
        "logs/eQTL_mapping/QQPlot.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/RealVsPermutatedQQPlot.R &> {log}
        """

rule GTExRemapLefflerSNPs:
    input:

