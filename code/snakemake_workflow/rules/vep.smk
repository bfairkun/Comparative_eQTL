rule vep_chimp_ref_panel:
    input:
    output:
        "test"
    log:
    conda:
        "../envs/vep.yaml"
    shell:
        """
        vep -i ../../PopulationSubstructure/ReferencePanelMerged.annotated.splits/22.vcf --format vcf -o Chimp.22.test.txt --stats_text --sf Chimp.22.test.stats.txt --cache --dir_cache cachedir/ --species pan_troglodytes --force --pick --coding_only
        """

rule vep_chimp_ref_panel_merge:

rule extract_test_snp_positions:
    """
    the vcf input is from converting vcf>plink>vcf which corroded the reference and
    variant alleles since plink re-encodes alleles based on allele frequency while
    the original vcf is based on reference/non-reference. Therefore, to get a vcf
    of test snps where the alleles are properly coded like the original vcf, I will
    use bcftools isec to intersect vcf files by position, writing out the alleles
    from the original vcf.
    """
    input:
        vcf_after_plink = "eQTL_mapping/FastQTL/ForAssociationTesting.vcf.gz",
        vcf_original = "PopulationSubstructure/ReferencePanelMerged.annotated.vcf.gz"
    output:
        vcf = "vep/vcf/Chimp.Ref_Plus_ThisStudy.vcf.gz",
        tbi = "vep/vcf/Chimp.Ref_Plus_ThisStudy.vcf.gz.tbi"
    shell:
        """
        bcftools isec -p vep/vcf -n=2 -O z -c all -w2 {input.vcf_after_plink} {input.vcf_original}
        mv vep/vcf/0001.vcf.gz {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule vep_chimp_tested_snps:

rule vep_chimp_tested_snps_merge:
