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
        grep 'ENSPTRG00000000104' eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt | awk -v OFS='\\t' -F'\\t' 'BEGIN{{print "MarkerName", "P-value"}} {{split($1,a,"."); print "chr"a[2]":"a[3], $5}}' | python2.7 /project/gilad/bjf79/software/locuszoom/bin/locuszoom --build PanTro5 --refgene ENSPTRG00000000104 --flank 250kb --metal - -p locuszoomtest  -v --source samples --pop MyPop --cache None
        """

# rule GetHumanFullEqtlResults:
#     output:
#     shell:
#         """
#         wget
#         """

# rule GetHumanLeadSnpsForMatchedWindowSize:

rule GetEqtlResultsForSharedPolymorphisms:
    input:
        eqtls = "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt",
        shared_polymorhisms = "../../data/LeflerSharedPolymorphisms.PanTro5.bed"
    output:
        shared = "../../output/SharedPolymorphisms.shared.chimpeqtls.txt",
        random = "../../output/SharedPolymorphisms.random.chimpeqtls.txt"
    shell:
        """
        sed 's/^chr//g' {input.shared_polymorhisms} | awk -F'\\t' '{{print "ID."$1"."$2}}' | grep -f - {input.eqtls} > {output.shared}
        awk 'NR>1' {input.eqtls} | shuf -n 1000 > {output.random}
        """
