rule GetReferencePopulationVcfChromList:
    input:
        config["ref"]["reference-population-vcf"],
    output:
        dynamic("MiscOutput/ReferencePopChromList/{chrom}"),
    shell:
        """
        (zcat {input} | awk '$1 !~/^#/ {{print $1}}' | sort | uniq | awk '{{print $1 > "MiscOutput/ReferencePopChromList/" $1 }}' ) &> logs/GetReferencePopulationVcfChromList.log
        """

rule SplitReferencePopulationVcfByChrom:
    input:
        vcf=config["ref"]["reference-population-vcf"],
        ChromList="MiscOutput/ReferencePopChromList/{chrom}"
    output:
        vcf = "PopulationSubstructure/ReferenceUnliftedByChrom/{chrom}.vcf.gz",
        tbi = "PopulationSubstructure/ReferenceUnliftedByChrom/{chrom}.vcf.gz.tbi"
    log:
        "logs/Misc/SplitReferencePopulationVcfByChrom/{chrom}.log"
    shell:
        """
        bcftools view -r {wildcards.chrom} -O z {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule LiftoverReferencePopulationVcf:
    """
    This step requires a lot of memory and will break for whole genome even if
    supplied the max amount of memory on broadwl. Therefore, it should be run
    one chromosome at a time. Hence the dynamic workflow.
    """
    input:
        vcf = "PopulationSubstructure/ReferenceUnliftedByChrom/{chrom}.vcf.gz"
    output:
        vcf = "PopulationSubstructure/ReferenceLiftedByChrom/{chrom}.vcf.gz",
        tbi = "PopulationSubstructure/ReferenceLiftedByChrom/{chrom}.vcf.gz.tbi"
    log:
        "logs/picard/LiftoverVcf/{chrom}.log"
    params:
        ref=config["ref"]["genome"],
        chain=config["ref"]["chain-file"]
    shell:
        """
        picard LiftoverVcf I={input.vcf} O={output.vcf} REJECT=/dev/null C={params.chain} R={params.ref} -Xmx45000m &> {log}
        tabix -f -p vcf {output.vcf} &>> {log}
        """

rule ConcatLiftedReferenceByChrom:
    input:
        vcf = dynamic("PopulationSubstructure/ReferenceLiftedByChrom/{chrom}.vcf.gz"),
        tbi = dynamic("PopulationSubstructure/ReferenceLiftedByChrom/{chrom}.vcf.gz.tbi")
    output:
        vcf = "PopulationSubstructure/ReferencePanel.catted.vcf.gz",
    log:
        "logs/Misc/ConcatLiftedReferenceByChrom.log"
    params:
    shell:
        """
        bcftools concat -a -O z {input.vcf} > {output.vcf} 2> {log}
        """

rule SortLiftedReference:
    input:
        "PopulationSubstructure/ReferencePanel.catted.vcf.gz"
    output:
        vcf = "PopulationSubstructure/ReferencePanel.vcf.gz",
    log:
        "logs/Misc/SortLiftedReference"
    params:
        "-m 10G -T /scratch/midway2/bjf79/TestResultsScratch/SortScratch/"
    shell:
        """
        bcftools sort {params} -O z {input} > {output.vcf} 2> {log}
        """


rule IndexLiftedReference:
    input:
        "PopulationSubstructure/ReferencePanel.vcf.gz"
    output:
        "PopulationSubstructure/ReferencePanel.vcf.gz.tbi"
    log:
        "logs/Misc/IndexLiftedReference"
    shell:
        """
        tabix -f -p vcf {input} &> {log}
        """

rule AddReferencePopulationToVCF:
    input:
        vcf = "filtered/all.vcf.gz",
        tbi = "filtered/all.vcf.gz.tbi",
        Refvcf = "PopulationSubstructure/ReferencePanel.vcf.gz",
        Reftbi = "PopulationSubstructure/ReferencePanel.vcf.gz.tbi",
    output:
        vcf = "PopulationSubstructure/ReferencePanelMerged.vcf",
        tbi = "PopulationSubstructure/ReferencePanelMerged.vcf.tbi",
    log:
        "logs/populationsubstructure/AddRefPopulationToVcf.log"
    shell:
        """
        bcftools merge {input.vcf} {input.Refvcf} -O z > {output.vcf} 2> {log}
        tabix -p vcf {output.vcf} &>> {log}
        """

rule MakePlinkBed:
    """
    bed is required for for admixture. tped is for REAP. geno 0 parameter to include snps with 0% missing genotypes
    """
    input:
        vcf = ancient("PopulationSubstructure/ReferencePanelMerged.vcf"),
    output:
        bed =  "PopulationSubstructure/plink/Merged.bed",
        fam = "PopulationSubstructure/plink/Merged.fam",
        bim = "PopulationSubstructure/plink/Merged.bim"
    log:
        "logs/plink/PopulationSubstructureMakeBed"
    params:
        config["PopulationSubstructure"]["AdmixturePlinkFilters"]
    shell:
        """
        plink --vcf {input.vcf} --vcf-half-call m --allow-extra-chr --out PopulationSubstructure/plink/Merged --geno 0 -recode12 {params} --make-bed &> {log}
        """

rule MakePlinkTped:
    """
    bed is required for for admixture. tped is for REAP. geno 0 parameter to include snps with 0% missing genotypes
    """
    input:
        vcf = ancient("PopulationSubstructure/ReferencePanelMerged.vcf"),
    output:
        tfam = "PopulationSubstructure/plink/Merged.tfam",
        tped =  "PopulationSubstructure/plink/Merged.tped",
    log:
        "logs/plink/PopulationSubstructureMakeTped"
    params:
        config["PopulationSubstructure"]["AdmixturePlinkFilters"]
    shell:
        """
        plink --vcf {input.vcf}  --vcf-half-call m --allow-extra-chr --out {output.base}  --geno 0 -recode12 {params} --transpose &> {log}
        """

rule FixChromsomeNamesForAdmixtureHack:
    """
    Admixture needs chromosome names as integers. 2A and 2B are not acceptable.
    This rule edits the bim file to make compatible. Also, copied bam and fam
    files to a matching prefix, since Admixture searches for matching prefixes.
    """
    input:
        bim = "PopulationSubstructure/plink/Merged.bim",
        bed = "PopulationSubstructure/plink/Merged.bed",
        fam = "PopulationSubstructure/plink/Merged.fam",
    output:
        bim = "PopulationSubstructure/plink/MergedForAdmixture.bim",
        bed = "PopulationSubstructure/plink/MergedForAdmixture.bed",
        fam = "PopulationSubstructure/plink/MergedForAdmixture.fam",
    log:
        "logs/Misc/FixChromsomeNamesForAdmixtureHack.log"
    shell:
        """
        cp {input.bed} {output.bed}
        cp {input.fam} {output.fam}
        sed 's/^2A\s/24\t/g; s/^2B\s/25\t/g' {input.bim} > {output.bim} 2> {log}
        """

rule Admixture:
    "admixture does not allow you to change default output filepaths, so output must be renamed manually"
    input:
        bed =  "PopulationSubstructure/plink/MergedForAdmixture.bed",
    output:
        P = "MergedForAdmixture.{K}.P",
        Q = "MergedForAdmixture.{K}.Q",
    log:
        "logs/Admixture/Admixture.{K}.log"
    params:
    shell:
        """
        admixture {input.bed} {wildcards.K} &> {log}
        """

rule ReformatAdmixture:
    input:
        P = "MergedForAdmixture.{K}.P",
        Q = "MergedForAdmixture.{K}.Q",
        fam = "PopulationSubstructure/plink/MergedForAdmixture.fam",
    output:
        P = "PopulationSubstructure/Admixture/MergedForAdmixture.{K}.P",
        Q = "PopulationSubstructure/Admixture/MergedForAdmixture.{K}.Q",
        Q_labelled = config["gitinclude_output"] + "PopulationStructure/Admixture/MergedForAdmixture.{K}.Q.labelled"
    log:
        "logs/Misc/ReformatAdmixture.{K}.log"
    shell:
        """
        mv {input.P} {output.P} 2>> {log}
        mv {input.Q} {output.Q} 2>> {log}
        paste {input.fam} {output.Q} > {output.Q_labelled} 2>> {log}
        """

rule PlotAdmixture:
    input:
        fam = "PopulationSubstructure/plink/MergedForAdmixture.fam",
        Q = "PopulationSubstructure/Admixture/MergedForAdmixture.{K}.Q"
    output:
        "PopulationSubstructure/Admixture/MergedForAdmixture.{K}.pdf"
    log:
        "logs/Misc/PlotAdmixture/{K}.log"
    shell:
        """
        Rscript scripts/plotAdmixture.R {input.Q} {input.fam} {output} &> {log}
        """

# rule REAP:
#     input:
#         P = "PopulationSubstructure/plink/Merged." + config["PopulationSubstructure"]["AdmixtureK"] + "P",
#         Q_labelled = "PopulationSubstructure/plink/Merged." + config["PopulationSubstructure"]["AdmixtureK"] + "Q.labelled",
#         tped =  "PopulationSubstructure/plink/Merged.tped",
#         tfam = "PopulationSubstructure/plink/Merged.tfam",
#     output:
#         "MyOutput"
#     log:
#         ""
#     params:
#         ""
#     shell:
#         """
#         ./REAP -g {input.tped} -p {input.tfam} -a {input.Q_labelled} -f {input.P} -r 1 -k {config[PopulationSubstructure][AdmixtureK]}
#         """

rule vcftools_relatedness2:
    input:
        "filtered/all.vcf.gz"
    output:
        "PopulationSubstructure/out.relatedness2"
    log:
        "logs/Misc/vcftools_relatedness2.log"
    params: "--maf 0.05 --remove-indv Pan_troglodytes_ThisStudy-4X0354 --remove-indv Pan_troglodytes_ThisStudy-462 --remove-indv Pan_troglodytes_ThisStudy-495 --remove-indv Pan_troglodytes_ThisStudy-4X0339 --remove-indv Pan_troglodytes_ThisStudy-554"
    shell:
        """
        vcftools --gzvcf {input} {params} --out PopulationSubstructure/out --relatedness2 &> {log}
        """

