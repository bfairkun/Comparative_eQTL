rule GetExonicRegions:
    """
    Create files to subset vcf with bcftools with vep (much faster than running
    vep on full vcf). Used CDS regions +/- 3bp to allow vep to search for
    alterations in stop codons which are 3bp outside of CDS boundaries.
    """
    input:
        gtf_chimp = "/project2/gilad/bjf79/genomes/Pan_tro_3.0_Ensembl/GeneAnnotation/Pan_troglodytes.Pan_tro_3.0.94.chr.gtf",
        gtf_human = "/project2/gilad/bjf79/genomes/GRCh38_Ensembl/Annotations/Homo_sapiens.GRCh38.94.chr.gtf"
    output:
        exonic_bed_chimp = "MiscOutput/ExonicRegions.chimp.bed",
        exonic_bed_human = "MiscOutput/ExonicRegions.human.bed"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$3=="CDS" {{ print $1, $4-4, $5+3 }}' {input.gtf_chimp} | bedtools sort -i - | bedtools merge -i - > {output.exonic_bed_chimp}
        awk -F'\\t' -v OFS='\\t' '$3=="CDS" {{ print "chr"$1, $4-4, $5+3 }}' {input.gtf_human} | bedtools sort -i - | bedtools merge -i - > {output.exonic_bed_human}
        """

rule vep:
    input:
        vcf = GetVcfForVep,
        tbi = GetTbiForVep,
        CDS_bed = "MiscOutput/ExonicRegions.{species}.bed"
    output:
        results = "vep/{species}/{IndvSubsample}/{chromosome}.txt",
    log:
        "logs/vep/{species}.{IndvSubsample}.{chromosome}.log"
    conda:
        "../envs/vep.yaml"
    params:
        ExtraBcftoolsParam = GetExtraParamsForBvctoolsBeforeVep,
        SpeciesDb = GetSpeciesDbKey,
        IntermediateCommand = GetCommandBetweenBcftoolsAndVep,
        chromprefix = GetChromPrefix
    wildcard_constraints:
        species="\w+",
        chromosome="\w+",
        IndvSubsample="\w+"
    shell:
        """
        bcftools view -q 0.1:minor --force-samples {params.ExtraBcftoolsParam}  -R <(cat {input.CDS_bed} | awk -F'\\t' '$1=="{params.chromprefix}{wildcards.chromosome}"') {input.vcf} {params.IntermediateCommand} | vep -i /dev/stdin --offline --format vcf -o {output.results} --no_stats --cache --dir_cache /project2/gilad/bjf79/genomes/vep_cache_dir/ --species {params.SpeciesDb} --force --pick --coding_only &> {log}
        """


rule CountPnPs:
    input:
        GetVepOutput
    output:
        # Pn = "vep/PnPs/{species}/{IndvSubsample}/Pn.txt",
        # Ps = "vep/PnPs/{species}/{IndvSubsample}/Ps.txt"
        # PnPs = "vep/PnPs/{species}/{IndvSubsample}/PnPs.txt"
        Catted = "vep/MergedResults/{species}/{IndvSubsample}/Merged.txt",
        PnPs = "../../output/NeutralityIndex/{species}/{IndvSubsample}/PnPs.txt"
    wildcard_constraints:
        species="\w+",
        IndvSubsample="\w+"
    shell:
        """
        cat {input} | awk '!/^#/' > {output.Catted}
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/PnPsMergeTables.R {output.Catted} {output.PnPs}
        """

        # cat {input} | awk -F'\\t' '!/^#/ {{ print $4, $7 }}' | awk '$2!~"synonymous" {{ print $1 }}' | sort | uniq -c | awk '{{ print $2, $1 }}' | sort > {output.Pn}
        # cat {input} | awk -F'\\t' '!/^#/ {{ print $4, $7 }}' | awk '$2~"synonymous" {{ print $1 }}' | sort | uniq -c | awk '{{ print $2, $1 }}' | sort > {output.Ps}
# rule vep_merge:

#     input:
#         expand("vep/{{species}}/{{IndvSubsample}}/{chrom}.txt", chrom=)

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
