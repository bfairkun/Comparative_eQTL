rule MakeSubsamples_ConstantSize:
    input: config["GTEx_eQTL_mapping"]["sample_list"]
    output:
        expand("GTEX_renalysis/Subsamples/SampleLists/{i}.txt", i=range(1, Num_subsamples + 1))
    params:
        n = config["GTEx_eQTL_mapping"]["subsample_size"]
    shell:
        """
        for i in $(seq 1 {Num_subsamples}); do
            shuf -n {params.n} {input} > GTEX_renalysis/Subsamples/SampleLists/${{i}}.txt
        done
        """

rule FastQTL_GTEx:
    input:
        samples = "GTEX_renalysis/Subsamples/SampleLists/{subsample}.txt",
        vcf = config["GTEx_eQTL_mapping"]["vcf"],
        vcf_tbi = config["GTEx_eQTL_mapping"]["vcf"] + ".tbi",
        phenotypes = config["GTEx_eQTL_mapping"]["phenotypes"],
        phenotypes_tbi = config["GTEx_eQTL_mapping"]["phenotypes"] + ".tbi",
        covariates = config["GTEx_eQTL_mapping"]["covariates"]
    output:
        permutation_chunk = "GTEX_renalysis/Subsamples/output_chunks/Subsample_{subsample}.Chunk_{chunk}.txt.gz"
    log:
        "logs/GTEx_reanalysis/FastQTL_GTEx/{subsample}.{chunk}.log"
    params:
        config["GTEx_eQTL_mapping"]["FastQTL_params"]
    shell:
        """
        fastQTL.1.165.linux --cov {input.covariates} --vcf {input.vcf} --include-samples {input.samples} --bed {input.phenotypes} --chunk {wildcards.chunk} {config[GTEx_eQTL_mapping][chunks]} --permute 1000 10000 {params} --out {output.permutation_chunk} &> {log}
        """

rule MergeFastQTL_GTEx:
    input:
        permutation_chunk = expand("GTEX_renalysis/Subsamples/output_chunks/Subsample_{{subsample}}.Chunk_{chunk}.txt.gz", chunk=range(1, 1+int(config["GTEx_eQTL_mapping"]["chunks"])))
    output:
        "GTEX_renalysis/Subsamples/output/Subsample_{subsample}.txt.gz"
    shell:
        """
        cat {input} > {output}
        """
