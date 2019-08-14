# rule MakeBestSnpBoxplots:

# rule AddChrToVcf:
#     input:
#         # vcf = "filtered/all.vcf.gz",
#         vcf = "eQTL_mapping/FastQTL/ForAssociationTesting.vcf.gz"
#     output:
#         ChromRename = "MiscOutput/ChromsomeList.rename.txt",
#         vcf = "MiscOutput/ForAssociationTesting.chromsrenamed.vcf.gz",
#         tabix = "MiscOutput/ForAssociationTesting.chromsrenamed.vcf.gz.tbi"
#     shell:
#         """
#         bcftools view -H {input.vcf} | awk -F'\\t' '{{print $1}}' |  uniq | awk '{{print $1, "chr"$1}}' > {output.ChromRename}
#         bcftools view {input.vcf} -i 'MAF[0]>0.1 & F_MISSING < 0.1' | bcftools annotate --rename-chrs {output.ChromRename} -O z {input.vcf} > {output.vcf}
#         tabix -p vcf {output.vcf}
#         """

rule AddChrToVcf:
    input:
        # vcf = "filtered/all.vcf.gz",
        vcf = "eQTL_mapping/FastQTL/ForAssociationTesting.vcf.gz"
    output:
        ChromRename = "MiscOutput/ChromsomeList.rename.txt",
        vcf = "MiscOutput/ForAssociationTesting.chromsrenamed.vcf.gz",
        tabix = "MiscOutput/ForAssociationTesting.chromsrenamed.vcf.gz.tbi"
    shell:
        """
        bcftools view -H {input.vcf} | awk -F'\\t' '{{print $1}}' |  uniq | awk '{{print $1, "chr"$1}}' > {output.ChromRename}
        bcftools annotate --rename-chrs {output.ChromRename} -O z {input.vcf} -x ID -I +'%CHROM:%POS' > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule PrepareChimpLocusZoomDb:
    input:
        gtf = config["ref"]["genomegtf"],
        snploc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        locuszoom_config = "../../data/locuszoom.config"
    output:
        genePred = "MiscOutput/Chimp.genepred",
        refflat = "MiscOutput/Chimp.refflat",
        locuszoomdb = config["LocusZoomDirectoryPath"] + "data/database/locuszoom_PanTro5.db",
        snploc = "MiscOutput/Chimp.snploc",
        snpset = "MiscOutput/Chimp.snpset",
        snptrans = "MiscOutput/Chimp.trans.txt",
        locuszoom_config = config["LocusZoomDirectoryPath"] + "conf/m2zfast.conf"
    log:
        "logs/LocusZoom/Makedatabase.log"
    shell:
        """
        gtfToGenePred -genePredExt {input.gtf} {output.genePred}

        # Convert to genepred by genes instead of transcripts; use longest transcript
        awk -F'\\t' -v OFS='\\t' '{{print $5-$4, $12, $0}}' {output.genePred} | sort -k2,2 -k1nr,1 | sort -u -k2,2 | awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "geneName", "name","chrom", "strand", "txStart","txEnd","cdsStart","cdsEnd", "exonCount","exonStarts","exonEnds"}}  {{print $2, $3, "chr"$4,$5,$6,$7,$8,$9,$10,$11, $12}}' > {output.refflat}

        # snp_pos file needs header with correct labels
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "snp", "chr", "pos"}} NR>1 {{print "chr"$2":"$3,"chr"$2, $3}}' {input.snploc} > {output.snploc}

        awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "rs_orig", "rs_current"}} NR>1 {{print $1, $1}}' {output.snploc} > {output.snptrans}

        # snp set file
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "snp", "snp_set"}} NR>1 {{print $1, "Illu1M"}}' {output.snploc} > {output.snpset}

        python2.7 {config[LocusZoomDirectoryPath]}bin/dbmeister.py  --db {output.locuszoomdb} --refflat {output.refflat} --snp_pos {output.snploc} --snp_set {output.snpset} --trans {output.snptrans} &> {log}

        mv {input.locuszoom_config} {output.locuszoom_config}
        """



rule PrepareHumanLocusZoomDb:
    input:
        gtf = config["ref"]["genomegtf"],
        snploc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        locuszoom_config = "../../data/locuszoom.config"
    output:
        genePred = "MiscOutput/Chimp.genepred",
        refflat = "MiscOutput/Chimp.refflat",
        locuszoomdb = config["LocusZoomDirectoryPath"] + "data/database/locuszoom_PanTro5.db",
        snploc = "MiscOutput/Chimp.snploc",
        snpset = "MiscOutput/Chimp.snpset",
        snptrans = "MiscOutput/Chimp.trans.txt",
        locuszoom_config = config["LocusZoomDirectoryPath"] + "conf/m2zfast.conf"
    log:
        "logs/LocusZoom/Makedatabase.log"
    shell:
        """
        gtfToGenePred -genePredExt {input.gtf} {output.genePred}

        # Convert to genepred by genes instead of transcripts; use longest transcript
        awk -F'\\t' -v OFS='\\t' '{{print $5-$4, $12, $0}}' {output.genePred} | sort -k2,2 -k1nr,1 | sort -u -k2,2 | awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "geneName", "name","chrom", "strand", "txStart","txEnd","cdsStart","cdsEnd", "exonCount","exonStarts","exonEnds"}}  {{print $2, $3, "chr"$4,$5,$6,$7,$8,$9,$10,$11, $12}}' > {output.refflat}

        # snp_pos file needs header with correct labels
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "snp", "chr", "pos"}} NR>1 {{print "chr"$2":"$3,"chr"$2, $3}}' {input.snploc} > {output.snploc}

        awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "rs_orig", "rs_current"}} NR>1 {{print $1, $1}}' {output.snploc} > {output.snptrans}

        # snp set file
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "snp", "snp_set"}} NR>1 {{print $1, "Illu1M"}}' {output.snploc} > {output.snpset}

        python2.7 {config[LocusZoomDirectoryPath]}bin/dbmeister.py  --db {output.locuszoomdb} --refflat {output.refflat} --snp_pos {output.snploc} --snp_set {output.snpset} --trans {output.snptrans} &> {log}

        mv {input.locuszoom_config} {output.locuszoom_config}
        """


rule PlinkFilesForLocusZoom:
    input:
        "eQTL_mapping/plink/ForAssociationTesting.bed",
        "eQTL_mapping/plink/ForAssociationTesting.fam"
    output:
        bed = config["LocusZoomDirectoryPath"] + "data/genotypes/MyPop/chr{chromosome}.bed",
        bim = config["LocusZoomDirectoryPath"] + "data/genotypes/MyPop/chr{chromosome}.bim"
    log:
        "logs/LocusZoom/MakeBedfiles/{chromosome}.log"
    shell:
        """
        plink --make-bed --bfile eQTL_mapping/plink/ForAssociationTesting --allow-extra-chr --chr {wildcards.chromosome} --out {config[LocusZoomDirectoryPath]}data/genotypes/MyPop/chr{wildcards.chromosome}
        cat {config[LocusZoomDirectoryPath]}data/genotypes/MyPop/chr{wildcards.chromosome}.bim | awk -F'\\t' -v OFS='\\t' '{{split($2,a,"."); print $1, "chr"$1":"a[3], $3,$4,$5,$6}}' > {config[LocusZoomDirectoryPath]}data/genotypes/MyPop/chr{wildcards.chromosome}.bim.temp
        mv {config[LocusZoomDirectoryPath]}data/genotypes/MyPop/chr{wildcards.chromosome}.bim.temp {output.bim}
        """

rule HumanPlinkFilesForLocusZoom:
    input:
        bim = "/project2/yangili1/GTEx/genotypes_bfile/GTEx_geno_{chromosome}.bim",
        bed = "/project2/yangili1/GTEx/genotypes_bfile/GTEx_geno_{chromosome}.bed",
        fam = "/project2/yangili1/GTEx/genotypes_bfile/GTEx_geno_{chromosome}.fam",
        plog = "/project2/yangili1/GTEx/genotypes_bfile/GTEx_geno_{chromosome}.log",
        nosex = "/project2/yangili1/GTEx/genotypes_bfile/GTEx_geno_{chromosome}.nosex"
    output:
        bim = config["LocusZoomDirectoryPath"] + "data/genotypes/gtex/chr{chromosome}.bim",
        bed = config["LocusZoomDirectoryPath"] + "data/genotypes/gtex/chr{chromosome}.bed",
        fam = config["LocusZoomDirectoryPath"] + "data/genotypes/gtex/chr{chromosome}.fam",
        plog = config["LocusZoomDirectoryPath"] + "data/genotypes/gtex/chr{chromosome}.log",
        nosex = config["LocusZoomDirectoryPath"] + "data/genotypes/gtex/chr{chromosome}.nosex",
    shell:
        """
        awk -F'\\t' -v OFS='\\t'  '{{ print "chr"$1, "chr"$1":"$4, $3, $4, $5,$6 }} ' {input.bim} > {output.bim}
        cp {input.bed} {output.bed}
        cp {input.fam} {output.fam}
        cp {input.plog} {output.plog}
        cp {input.nosex} {output.nosex}
        """
 
# rule LocusZoom:
#     input:
#         locuszoomdb = "MiscOutput/Chimp.locuszoom.db",
#         expand(config["LocusZoomDirectoryPath"] + "data/genotypes/MyPop/chr{chromosome}", chromosome=[str(i) for i in range(3,23)] + ['1', '2A', '2B'])
#     output:
#         "MiscOutput/LocusZoomPlot"
#     shell:
#         """
#         grep 'ENSPTRG00000000104' eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt | awk -v OFS='\\t' -F'\\t' 'BEGIN{{print "MarkerName", "P-value"}} {{split($1,a,"."); print "chr"a[2]":"a[3], $5}}' | python2.7 /project/gilad/bjf79/software/locuszoom/bin/locuszoom --build PanTro5 --refgene ENSPTRG00000000104 --flank 250kb --metal - -p locuszoomtest --ld-vcf {input.vcf} -v  --ignore-vcf-filter --cache None
#         """

rule LocusZoom:
    input:
        locuszoomdb = config["LocusZoomDirectoryPath"] + "data/database/locuszoom_PanTro5.db",
        plink_files_for_ld = expand(config["LocusZoomDirectoryPath"] + "data/genotypes/MyPop/chr{chromosome}.bed", chromosome=[str(i) for i in range(3,23)] + ['1', '2A', '2B'])
    output:
        "MiscOutput/LocusZoomPlot"
    shell:
        """
        grep 'ENSPTRG00000016090' eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt | awk -v OFS='\\t' -F'\\t' 'BEGIN{{print "MarkerName", "P-value"}} {{split($1,a,"."); print "chr"a[2]":"a[3], $5}}' | python2.7 /project/gilad/bjf79/software/locuszoom/bin/locuszoom --build PanTro5 --refgene ENSPTRG00000016090 --flank 250kb --metal - -p locuszoomtest  -v --source SOURCE --pop MyPop --cache None --bed-tracks ../../data/LeflerSharedPolymorphisms.PanTro5.bed
        """

# Plot two reference snps: top snp, and top human shared snp
# grep 'ENSPTRG00000016090' eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt | awk -v OFS='\t' -F'\t' 'BEGIN{print "MarkerName", "P-value"} {split($1,"."); print "chr"a[2]":"a[3], $5}' | python2.7 /project/gilad/bjf79/software/locuszoom/bin/locuszoom --build PanTro5 --refsnp chr4:73517684 --flank 250kb --metal - -p locuszoomtest  -v --source SOURCE --pop MyPop --cache None --add-refsnps "chr4:73468307" --bed ../../data/LeflerSharedPolymorphisms.PanTro5.bed)}'

rule LocusZoomHuman:
    input:
        plink_files_for_ld = expand(config["LocusZoomDirectoryPath"] + "data/genotypes/gtex/chr{chromosome}.bed", chromosome=[str(i) for i in range(1,23)]),
        locuszoom_config = config["LocusZoomDirectoryPath"] + "conf/m2zfast.conf"
    output:
        "MiscOutput/HumanLocusZoomPlot"
    shell:
        """
        cat ../../data/HumanGeneSnpStats.ENSG00000163453.hg38.txt | awk -v OFS='\\t' -F'\\t' 'BEGIN{{print "MarkerName", "P-value"}} {{split($2,"_"); print "chr"a[1]":"a[2], $7}}' | python2.7 /project/gilad/bjf79/software/locusZoomWithLD/locuszoom/bin/locuszoom --build hg38 --flank 250kb --refsnp chr4:57913756 --metal - -v --source 1000G_Nov2014 --pop EUR --cache None --plotonly --no-date -p MiscOutput/MyGene'
        """

# rule GetHumanFullEqtlResults:
#     output:
#     shell:
#         """
#         wget
#         """

# rule GetHumanLeadSnpsForMatchedWindowSize:

# rule GetEqtlResultsForSharedPolymorphisms:
#     input:
#         eqtls = "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt",
#         shared_polymorhisms = "../../data/LeflerSharedPolymorphisms.PanTro5.bed"
#     output:
#         shared = "../../output/SharedPolymorphisms.shared.chimpeqtls.txt",
#         random = "../../output/SharedPolymorphisms.random.chimpeqtls.txt"
#     shell:
#         """
#         sed 's/^chr//g' {input.shared_polymorhisms} | awk -F'\\t' '{{print "ID."$1"."$2}}' | grep -f - {input.eqtls} > {output.shared}
#         awk 'NR>1' {input.eqtls} | shuf -n 1000 > {output.random}
#         """

rule FixHumanVcfHeaders_YRI:
    input:
        lambda wildcards: "/project2/yangili1/ankeetashah/noisysQTL/EU_LCL/vcf_2/chr{chromsome}.keep.recode.two.vcf.gz".format(chromsome = wildcards.chromosome)
    output:
        vcf = "MiscOutput/HumanVcfs/ByChrom/YRI.{chromosome}.vcf.gz",
        tbi = "MiscOutput/HumanVcfs/ByChrom/YRI.{chromosome}.vcf.gz.tbi"
    log:
        "logs/Misc/FixHumanVcfHeaders_YRI/{chromosome}.log"
    shell:
        """
        zcat {input} | awk 'NR==2 {{printf("##FORMAT=<ID=BD,Number=1,Type=String,Description=\\"?\\">\\n##FORMAT=<ID=PP,Number=1,Type=String,Description=\\"?\\">\\n");}} {{print}}' | bcftools annotate -h ../../data/hg19.chromsizes.vcf.headerlines -O z - > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule ConcatHumanVcfHeaders_YRI:
    input:
        vcf = expand("MiscOutput/HumanVcfs/ByChrom/YRI.{chromosome}.vcf.gz", chromosome=[str(i) for i in range(1,23)])
    output:
        vcf = "MiscOutput/HumanVcfs/YRI.vcf.gz",
        tbi = "MiscOutput/HumanVcfs/YRI.vcf.gz.tbi"
    log:
        "logs/Misc/ConcatHumanVcfHeaders_YRI.log"
    shell:
        """
        bcftools concat -a -O z {input.vcf} > {output.vcf} 2> {log}
        tabix -p vcf {output.vcf}
        """


rule FixHumanVcfHeaders_GEUVADIS:
    input:
        lambda wildcards: "/project2/yangili1/ankeetashah/noisysQTL/EU_LCL/vcf_2/chr{chromsome}.keep.recode.two.vcf.gz".format(chromsome = wildcards.chromosome)
    output:
        vcf = "MiscOutput/HumanVcfs/ByChrom/GEUVADIS.{chromosome}.vcf.gz",
        tbi = "MiscOutput/HumanVcfs/ByChrom/GEUVADIS.{chromosome}.vcf.gz.tbi"
    log:
        "logs/Misc/FixHumanVcfHeaders_GEUVADIS/{chromosome}.log"
    shell:
        """
        zcat {input} | awk 'NR==2 {{printf("##FORMAT=<ID=BD,Number=1,Type=String,Description=\\"?\\">\\n##FORMAT=<ID=PP,Number=1,Type=String,Description=\\"?\\">\\n");}} {{print}}' | bcftools annotate -h ../../data/hg19.chromsizes.vcf.headerlines -O z - > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule ConcatHumanVcfHeaders_GEUVADIS:
    input:
        vcf = expand("MiscOutput/HumanVcfs/ByChrom/GEUVADIS.{chromosome}.vcf.gz", chromosome=[str(i) for i in range(1,23)])
    output:
        vcf = "MiscOutput/HumanVcfs/GEUVADIS.vcf.gz",
        tbi = "MiscOutput/HumanVcfs/GEUVADIS.vcf.gz.tbi"
    log:
        "logs/Misc/ConcatHumanVcfHeaders_GEUVADIS.log"
    shell:
        """
        bcftools concat -a -O z {input.vcf} > {output.vcf} 2> {log}
        tabix -p vcf {output.vcf}
        """

rule LDDecayPlots:
    input:
        H_geuvadis = "MiscOutput/HumanVcfs/GEUVADIS.vcf.gz",
        H_YRI = "MiscOutput/HumanVcfs/YRI.vcf.gz",
        C = "eQTL_mapping/FastQTL/ForAssociationTesting.vcf.gz"
    output:
        H_geuvadis = "MiscOutput/LDDecay/H_geuvadis.stat.gz",
        H_YRI = "MiscOutput/LDDecay/H_YRI.stat.gz",
        C = "MiscOutput/LDDecay/C.stat.gz"
    shell:
        """
        /project/gilad/bjf79/software/PopLDdecay/bin/PopLDdecay -InVCF {input.C} -OutStat {output.C}
        /project/gilad/bjf79/software/PopLDdecay/bin/PopLDdecay -InVCF {input.H_YRI} -OutStat {output.H_YRI}
        /project/gilad/bjf79/software/PopLDdecay/bin/PopLDdecay -InVCF {input.H_geuvadis} -OutStat {output.H_geuvadis}
        """

rule GetAllSnpGenePairsForLeflerSharedPolymorphisms:
    input:
        SharedPolymorphisms = "../../data/LeflerSharedPolymorphisms.PanTro5.bed",
        Vcf = "MiscOutput/ForAssociationTesting.chromsrenamed.vcf.gz",
        eqtl_results = "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt",
        chromsizes = "../../data/PanTro5.chrome.sizes"
    output:
        Shared = "../../output/SharedPolymorphisms.shared.chimpeqtls.txt",
        Random = "../../output/SharedPolymorphisms.random.chimpeqtls.txt"
    shell:
        """
        bedtools slop -i {input.SharedPolymorphisms} -r 1 -l 1 -g {input.chromsizes} | bcftools view -H -R - {input.Vcf} | awk -F'\\t' '{{print $1"."$2}}' | sed -e 's/^chr/ID\./' | grep -w -f - {input.eqtl_results} > {output.Shared}
        awk 'NR>1' {input.eqtl_results} | shuf -n 1000 > {output.Random}
        """

rule DownloadGTExSummaryStatsAllTissues:
    output:
        expand()
    log:
        "logs/DownloadGTExSummaryStats.log"
    shell:
        """
        wget ...*
        """


rule MergeGTExSummaryStatsAllTissues:
    input:
        expand()
    log:
        "logs/CombineGTExSummaryStats.log"
    output:
        ""
    shell:
        """

        """
