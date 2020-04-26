# fastQTL --vcf genotypes.vcf.gz --bed phenotypes.bed.gz --region 22:17000000-18000000 --permute 1000 --out permutations.default.txt.gzfastQTL --vcf genotypes.vcf.gz --bed phenotypes.bed.gz --region 22:17000000-18000000 --permute 1000 --out permutations.default.txt.gz

rule DownloadGtexData:
    output:
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Uterus.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Thyroid.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Hippocampus.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Esophagus_Muscularis.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Minor_Salivary_Gland.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Cerebellar_Hemisphere.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Colon_Transverse.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Kidney_Cortex.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Lung.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Caudate_basal_ganglia.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Ovary.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Artery_Aorta.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Pancreas.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Cerebellum.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Pituitary.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Prostate.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Heart_Left_Ventricle.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Artery_Tibial.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Cells_EBV-transformed_lymphocytes.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Adrenal_Gland.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Testis.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Uterus.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Stomach.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Adipose_Subcutaneous.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Cortex.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Colon_Sigmoid.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Nerve_Tibial.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Frontal_Cortex_BA9.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Esophagus_Muscularis.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Spinal_cord_cervical_c-1.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Skin_Not_Sun_Exposed_Suprapubic.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Hypothalamus.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Hippocampus.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Artery_Tibial.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Minor_Salivary_Gland.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Esophagus_Gastroesophageal_Junction.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Skin_Not_Sun_Exposed_Suprapubic.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Liver.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Anterior_cingulate_cortex_BA24.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Putamen_basal_ganglia.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Cerebellar_Hemisphere.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Heart_Left_Ventricle.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Muscle_Skeletal.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Anterior_cingulate_cortex_BA24.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Stomach.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Ovary.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Pituitary.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Cortex.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Spinal_cord_cervical_c-1.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Adrenal_Gland.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Colon_Transverse.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Vagina.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Kidney_Cortex.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Heart_Atrial_Appendage.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Spleen.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Cells_Cultured_fibroblasts.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Frontal_Cortex_BA9.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Nucleus_accumbens_basal_ganglia.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Heart_Atrial_Appendage.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Esophagus_Mucosa.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Artery_Coronary.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Liver.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Cells_Cultured_fibroblasts.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Esophagus_Gastroesophageal_Junction.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Lung.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Thyroid.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Cerebellum.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Substantia_nigra.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Nucleus_accumbens_basal_ganglia.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Amygdala.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Putamen_basal_ganglia.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Muscle_Skeletal.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Skin_Sun_Exposed_Lower_leg.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Artery_Aorta.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Amygdala.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Nerve_Tibial.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Colon_Sigmoid.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Prostate.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Adipose_Visceral_Omentum.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Testis.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Esophagus_Mucosa.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Caudate_basal_ganglia.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Spleen.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Adipose_Visceral_Omentum.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Pancreas.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Substantia_nigra.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Small_Intestine_Terminal_Ileum.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Skin_Sun_Exposed_Lower_leg.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Breast_Mammary_Tissue.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Cells_EBV-transformed_lymphocytes.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Adipose_Subcutaneous.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Breast_Mammary_Tissue.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Vagina.v8.normalized_expression.bed.gz",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Brain_Hypothalamus.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Small_Intestine_Terminal_Ileum.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Artery_Coronary.v8.normalized_expression.bed.gz.tbi",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Esophagus_Gastroesophageal_Junction.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Stomach.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Minor_Salivary_Gland.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Thyroid.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Heart_Atrial_Appendage.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Adipose_Subcutaneous.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Kidney_Cortex.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Colon_Sigmoid.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Breast_Mammary_Tissue.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Nerve_Tibial.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Nucleus_accumbens_basal_ganglia.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Liver.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Testis.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Adrenal_Gland.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Esophagus_Mucosa.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Cerebellar_Hemisphere.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Lung.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Artery_Aorta.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Substantia_nigra.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Artery_Coronary.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Colon_Transverse.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Artery_Tibial.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Hippocampus.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Cerebellum.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Cells_Cultured_fibroblasts.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Amygdala.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Skin_Not_Sun_Exposed_Suprapubic.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Hypothalamus.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Vagina.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Prostate.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Skin_Sun_Exposed_Lower_leg.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Pancreas.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Anterior_cingulate_cortex_BA24.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Putamen_basal_ganglia.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Frontal_Cortex_BA9.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Adipose_Visceral_Omentum.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Uterus.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Ovary.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Spleen.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Spinal_cord_cervical_c-1.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Heart_Left_Ventricle.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Caudate_basal_ganglia.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Brain_Cortex.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Small_Intestine_Terminal_Ileum.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Cells_EBV-transformed_lymphocytes.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Pituitary.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Esophagus_Muscularis.v8.covariates.txt",
        "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Muscle_Skeletal.v8.covariates.txt",
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

