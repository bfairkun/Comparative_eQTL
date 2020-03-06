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

rule GetTestSNPsForEgenes:
    input:
        eqtls = "../../output/ChimpEgenes.eigenMT.txt.gz",
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        genelocs = "../../output/Genes.bed",
        fai = config["ref"]["genome"] + ".fai"
    output:
        "../../output/ChimpEgenes.testsnps.txt.gz"
    params:
        window = 250000
    shell:
        """
        zcat {input.eqtls} | awk -F'\\t' 'NR>1 && $6 <=0.1 {{print $2}}' | grep -w -f - {input.genelocs} | bedtools slop -b {params.window} -i - -g {input.fai} | bedtools intersect -a <(awk -F'\\t' -v OFS='\\t' 'NR>1 {{print $2, $3, $3+1}}' {input.snps}) -b - -wo | awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "chr", "snp.pos", "gene", "strand", "TSS.dist", "TES.dist"}} $9=="+" {{print $1,$2,$7,$9, $2-($5+250000), $2-($6-250000)}} $9=="-" {{print $1,$2,$7,$9,($6-250000)-$2, ($5+250000)-$2}}' | gzip - > {output}
        """

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

rule index_YRI_fromKenneth:
    input: "/project2/gilad/kenneth/YRI/VCFs/genotypes/ALL.chr{MyChrom}_GRCh38.genotypes.20170504.vcf.gz"
    output: "/project2/gilad/kenneth/YRI/VCFs/genotypes/ALL.chr{MyChrom}_GRCh38.genotypes.20170504.vcf.gz.tbi"
    shell: "tabix -p vcf {input}"

rule ConcatHumanVcf_YRI_fromKenneth:
    input:
        vcf = expand("/project2/gilad/kenneth/YRI/VCFs/genotypes/ALL.chr{MyChrom}_GRCh38.genotypes.20170504.vcf.gz", MyChrom=[str(i) for i in range(1,23)] + ['X']),
        tbi =  expand("/project2/gilad/kenneth/YRI/VCFs/genotypes/ALL.chr{MyChrom}_GRCh38.genotypes.20170504.vcf.gz.tbi", MyChrom=[str(i) for i in range(1,23)] + ['X'])
    output:
        vcf = "MiscOutput/HumanVcfs/YRI.hg38.vcf.gz",
        tbi = "MiscOutput/HumanVcfs/YRI.hg38.vcf.gz.tbi"
    log:
        "logs/Misc/ConcatHumanVcf_YRI_fromKenneth.log"
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

rule GetPlinkMAF:
    input:
        bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
    output:
        frq = "MiscOutput/ForAssociationTesting.frq",
        zipped = "MiscOutput/allele_encodings.txt.gz"
    shell:
        """
        plink --freq --bfile eQTL_mapping/plink/ForAssociationTesting --allow-extra-chr --out MiscOutput/ForAssociationTesting
        cat {output.frq} | awk -v OFS='\\t' '{{ print $2,$3,$4,$5 }}' | gzip - > {output.zipped}
        """

rule GetSnpsMatchedToLeflerSharedPolymorphisms:
    input:
        SharedPolymorphisms = "../../data/LeflerSharedPolymorphisms.PanTro5.bed", 
        Vcf = "eQTL_mapping/FastQTL/ForAssociationTesting.vcf.gz",
        tbi = "eQTL_mapping/FastQTL/ForAssociationTesting.vcf.gz.tbi",
        chromsizes = "../../data/PanTro5.chrome.sizes",
        bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
    output:
        snps = "MiscOutput/SpeciesSharedSNPs.snps",
        shared_snp_loc = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedSnps.PanTro5.snploc",
        LD = "MiscOutput/SpeciesSharedSNPs.ld",
        snp_loc = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedMatchedSnps.PanTro5.snploc"
    shell:
        """
        bedtools slop -i {input.SharedPolymorphisms} -r 1 -l 1 -g {input.chromsizes} | sed 's/^chr//' | bcftools view -H -R - {input.Vcf} | awk -F'\\t' '{{ print $3 }}' | tee {output.snps} | awk -v OFS='\\t' 'BEGIN {{ print "snp", "chrom", "pos" }} {{ split($1,a,"."); print $1, a[2],a[3] }}' > {output.shared_snp_loc}
        plink --r2 with-freqs --bfile eQTL_mapping/plink/ForAssociationTesting --ld-snp-list {output.snps} --ld-window-kb 100 --ld-window-r2 0 --ld-window 99999 --allow-extra-chr --out MiscOutput/SpeciesSharedSNPs
        Rscript scripts/GetLeflerMatchedSnps.R {output.LD} {output.snp_loc}
        """

rule GetSpeciesSharedSNPGenotypesForMatrixEqtl:
    input:
        shared_snps  = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedSnps.PanTro5.snploc",
        shared_snps_matched = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedMatchedSnps.PanTro5.snploc",
        genotypes = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps"
    output:
        shared_snps_matched = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedMatchedSnps.PanTro5.snps",
        shared_snps = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedSnps.PanTro5.snps",
    shell:
        """
        cat <(head -1 {input.genotypes}) <(awk -F'\\t' 'NR>1 {{ print $1 }}' {input.shared_snps} | grep -w -f - {input.genotypes}) > {output.shared_snps}
        cat <(head -1 {input.genotypes}) <(awk -F'\\t' 'NR>1 {{ print $1 }}' {input.shared_snps_matched} | grep -w -f - {input.genotypes}) > {output.shared_snps_matched}
        """

rule MatrixEQTL_BestModelFromConfigFullResults_shared:
    """
    Matrix eQTL with full output for every snp-gene pair.
    """
    input:
        snps = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedSnps.PanTro5.snps",
        snp_locs = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedSnps.PanTro5.snploc",
        phenotypes = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt",
        gene_loc = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt",
        covariates = config["eQTL_mapping"]["CovariatesForFullOutput"],
        GRM = CovarianceMatrix,
    output:
        results = "eQTL_mapping/SharedPolymorphisms/MatrixEQTL.shared.results.txt",
        fig = "eQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/images/Results.png",
    log:
        "logs/eQTL_mapping/MatrixEQTL/ConfigCovariateModel.log"
    shell:
        """
        Rscript scripts/MatrixEqtl_Cis.AllPvals.R {input.snps} {input.snp_locs} {input.phenotypes} {input.gene_loc} {input.covariates} {input.GRM} {output.results} {output.fig} 250000 &> {log}
        """

rule wget_Gtex_tissue_MedianExpression:
    output:
        "../../data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz"
    shell:
        "wget -O {output} https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz"

rule CalculateTissueSpecificity:
    input:
        "../../data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz"
    output:
        TissueMatrixCondensedBrain = "Misc/CalculateTissueSpecificity/Matrix.CondensedBrain.tsv",
        tau = "../../output/TissueSpecificity/tau.txt",
        gini = "../../output/TissueSpecificity/gini.txt",
        tau_log = "../../output/TissueSpecificity/tau.log.txt",
        gini_log = "../../output/TissueSpecificity/gini.log.txt",
    shell:
        """
        Rscript scripts/SubsetGtexTissueMatrix.R {input} {output.TissueMatrixCondensedBrain}
        tspex {output.TissueMatrixCondensedBrain} {output.tau} tau
        tspex --log {output.TissueMatrixCondensedBrain} {output.tau_log} tau
        tspex {output.TissueMatrixCondensedBrain} {output.gini} gini
        tspex --log {output.TissueMatrixCondensedBrain} {output.gini_log} gini
        """

rule CalculateHeartCellTypeSpecificity:
    input:
        "../../data/CIBERSORT.Reference_sample_file.tsv"
    output:
        tau = "../../output/TissueSpecificity/CellTypeSpecificity.TabulaMurisHeart.tau.log.txt",
        z = "../../output/TissueSpecificity/CellTypeSpecificity.TabulaMurisHeart.z.log.txt"
    shell:
        """
        tspex --log {input} {output.tau} tau
        tspex --log {input} {output.z} zscore
        """

rule ReorganizeAndLiftLeflerSnps:
    input:
        TableS = "../../data/Lefler_TableS4.txt",
        hg19Tohg38Chain = "/project2/gilad/bjf79/genomes/LiftoverChains/hg19ToHg38.over.chain.gz",
        hg38ToPanTro5Chain = "/project2/gilad/bjf79/genomes/LiftoverChains/hg38ToPanTro5.over.chain.gz",
        GtexVcf = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz",
        covariates = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Heart_Left_Ventricle.v8.covariates.txt",
        ChimpVcf = "PopulationSubstructure/ReferencePanelMerged.annotated.vcf.gz",
    output:
        SharedSnpsHg38 = "MiscOutput/Lefler.hg38.bed",
        SharedSnpsHg38Named = "MiscOutput/Lefler.hg38.list.txt",
        SharedSnpsHg38NamedFreq = "MiscOutput/Lefler.hg38.frq",
        SharedSnpsPanTro5 = "MiscOutput/Lefler.PanTro5.bed",
        SharedSnpsPanTro5Named = "MiscOutput/Lefler.PanTro5.list.txt",
        SharedSnpsPanTro5NamedFreq = "MiscOutput/Lefler.PanTro5.frq",
        ClosestMatchingSnps = "MiscOutput/Lefler.PanTro5.rsID.chimpSNP.tab",
        CombinedOutput = "../../output/LeflerTestedSnps.tsv"
    shell:
        """
        #Liftover Lefler data to hg38
        cat {input.TableS} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{print "chr"$5,$6,$6+1, $9"."$13"."$14}}' | liftOver /dev/stdin {input.hg19Tohg38Chain} /dev/stdout /dev/null > {output.SharedSnpsHg38}
        #Grab SNPs with human SNP name
        bcftools view -S <(cat {input.covariates} | head -1 |transpose | awk 'NR>1') -q 0.01:minor -H -R <(cat {output.SharedSnpsHg38}) {input.GtexVcf} | awk -F'\\t' '{{ print $1,$2,$3 }}' > {output.SharedSnpsHg38Named}
        #Get allele frequency in this cohort
        bcftools view -S <(cat {input.covariates} | head -1 |transpose | awk 'NR>1') -q 0.01:minor -R <(cat {output.SharedSnpsHg38}) {input.GtexVcf} | vcftools --vcf - --freq --out MiscOutput/Lefler.hg38
        #Liftover to PanTro5
        liftOver {output.SharedSnpsHg38} {input.hg38ToPanTro5Chain} /dev/stdout /dev/null | sed 's/^chr//' > {output.SharedSnpsPanTro5}
        #Get SNPs with chimp SNP name
        bcftools view -S <(bcftools query -l {input.ChimpVcf} | grep "ThisStudy" | grep -v "MD_And") -q 0.1:minor -H -R <(cat {output.SharedSnpsPanTro5} | awk -F'\\t' -v OFS='\\t' '{{ print $1, $2-1, $3+1, $4 }}') {input.ChimpVcf} | awk -F'\\t' '{{ print $1,$2,$3 }}' > {output.SharedSnpsPanTro5Named}
        #Get AF in the chimp cohort
        bcftools view -S <(bcftools query -l {input.ChimpVcf} | grep "ThisStudy" | grep -v "MD_And") -q 0.1:minor -R <(cat {output.SharedSnpsPanTro5} | awk -F'\\t' -v OFS='\\t' '{{ print $1, $2-1, $3+1, $4 }}') {input.ChimpVcf} | vcftools --vcf - --freq --out MiscOutput/Lefler.PanTro5
        #Create translation table of snp positions where the snp position is off by 1.
        bedtools closest -a <(bedtools sort -i {output.SharedSnpsPanTro5}) -b <(awk -v OFS='\\t' '{{print $1, $2, $2+1, $3}}' {output.SharedSnpsHg38Named} | bedtools sort -i -) -d | awk -F'\\t' -v OFS='\\t' '{{split($4,H,"."); split($8,C,"."); print $4,$8, $NF, H[2]H[3], H[3]H[2], C[4]C[5],C[5]C[4]}}' | awk -F'\\t' -v OFS='\\t' '$3<3' | awk '$4==$6 || $4==$7' > {output.ClosestMatchingSnps}
        
        #Rscript to combine relevant info
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/ReorganizeLeflerSharedSnps.R
        """

rule PrepareLeflerSnpsForMatrixEQTL:
    input:
        LeflerSnps = "../../data/LeflerSharedPolymorphisms.hg38.bed",
        GtexVcf = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz",
        GtexVcfTbi = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz.tbi",
        covariates = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Heart_Left_Ventricle.v8.covariates.txt",
        expression = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Heart_Left_Ventricle.v8.normalized_expression.bed.gz"
    output:
        GenotypeRaw = "GTEX_renalysis/LeflerSnps/Genotypes.012",
        IndvRaw = "GTEX_renalysis/LeflerSnps/Genotypes.012.indv",
        PosRaw = "GTEX_renalysis/LeflerSnps/Genotypes.012.pos",
        GenotypeFinished = "GTEX_renalysis/LeflerSnps/MatrixEQTL/Genotypes.txt",
        SnpPos = "GTEX_renalysis/LeflerSnps/MatrixEQTL/snploc.txt",
        expression = "GTEX_renalysis/LeflerSnps/MatrixEQTL/expression.txt",
        geneloc = "GTEX_renalysis/LeflerSnps/MatrixEQTL/geneloc.txt"
    shell:
        """
        bcftools view -S <(cat {input.covariates} | head -1 |transpose | awk 'NR>1') -q 0.01:minor -R <(cat {input.LeflerSnps}) {input.GtexVcf} | vcftools --vcf - --012 --out GTEX_renalysis/LeflerSnps/Genotypes
        paste -d '\\t' {output.IndvRaw} {output.GenotypeRaw} | transpose | awk -F'\\t' -v OFS='\\t' 'NR==1 {{ print "id", $0 }} NR>2 {{ print "snp_"NR-2, $0 }}' > {output.GenotypeFinished}
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "snp", "chr", "pos" }} {{ print "snp_"NR, $1, $2 }}' {output.PosRaw} > {output.SnpPos}
        zcat {input.expression} | awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "geneid", "chr", "s1", "s2"}} NR>1 {{print $4, $1, $2, $3}}' > {output.geneloc}
        zcat {input.expression} | transpose | awk 'NR>3' | transpose > {output.expression}
        """

rule GetHumanSnpsMatchedToLeflerSharedPolymorphisms:
    input:
        LeflerSnps = "../../data/LeflerSharedPolymorphisms.hg38.bed",
        GtexVcf = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz",
        GtexVcfTbi = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz.tbi",
        covariates = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Heart_Left_Ventricle.v8.covariates.txt",
    output:
        bed = "GTEX_renalysis/LeflerSnps/ControlMatrixEQTL/SnpRegions.bed",
        snps = "GTEX_renalysis/LeflerSnps/ControlMatrixEQTL/SharedSNPs.txt",
        LD = "GTEX_renalysis/LeflerSnps/ControlMatrixEQTL/SharedSNPs.ld",
        snp_loc_raw = "GTEX_renalysis/LeflerSnps/ControlMatrixEQTL/ControlSnps.loc"
    shell:
        """
        bcftools view -S <(cat {input.covariates} | head -1 |transpose | awk 'NR>1') -q 0.01:minor -H -R <(cat {input.LeflerSnps}) {input.GtexVcf} | awk -F'\\t' '{{ print $3 }}' > {output.snps}
        bcftools view -S <(cat {input.covariates} | head -1 |transpose | awk 'NR>1') -q 0.01:minor -R <(cat {input.LeflerSnps} | awk -F'\\t' -v OFS='\\t' '{{ print $1, $2-100000, $3+100000 }}' | bedtools merge -i - ) {input.GtexVcf} | bcftools norm -d all | plink --vcf /dev/stdin --vcf-half-call m --make-bed --out GTEX_renalysis/LeflerSnps/ControlMatrixEQTL/SnpRegions
        plink --r2 with-freqs --bfile GTEX_renalysis/LeflerSnps/ControlMatrixEQTL/SnpRegions --ld-snp-list {output.snps} --ld-window-kb 100 --ld-window-r2 0 --ld-window 99999 --allow-extra-chr --out GTEX_renalysis/LeflerSnps/ControlMatrixEQTL/SharedSNPs
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/GetLeflerMatchedSnps.R {output.LD} {output.snp_loc}
        """

rule PrepareLeflerHumanMatchedSnpsForMatrixEQTL:
    input:
        MatchedSnps = "GTEX_renalysis/LeflerSnps/ControlMatrixEQTL/ControlSnps.loc",
        GtexVcf = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz",
        GtexVcfTbi = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz.tbi",
        covariates = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Heart_Left_Ventricle.v8.covariates.txt",
        expression = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_expression_matrices/Heart_Left_Ventricle.v8.normalized_expression.bed.gz"
    output:
        GenotypeRaw = "GTEX_renalysis/LeflerSnps/Genotypes.Control.012",
        IndvRaw = "GTEX_renalysis/LeflerSnps/Genotypes.Control.012.indv",
        PosRaw = "GTEX_renalysis/LeflerSnps/Genotypes.Control.012.pos",
        GenotypeFinished = "GTEX_renalysis/LeflerSnps/MatrixEQTL/Genotypes.control.txt",
        SnpPos = "GTEX_renalysis/LeflerSnps/MatrixEQTL/snploc.control.txt",
    shell:
        """
        bcftools view -S <(cat {input.covariates} | head -1 |transpose | awk 'NR>1') -q 0.01:minor -R <(cat {input.MatchedSnps} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{ print "chr"$2, $3, $3+1 }}') {input.GtexVcf} | vcftools --vcf - --012 --out GTEX_renalysis/LeflerSnps/Genotypes.Control
        paste -d '\\t' {output.IndvRaw} {output.GenotypeRaw} | transpose | awk -F'\\t' -v OFS='\\t' 'NR==1 {{ print "id", $0 }} NR>2 {{ print "snp_"NR-2, $0 }}' > {output.GenotypeFinished}
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "snp", "chr", "pos" }} {{ print "snp_"NR, $1, $2 }}' {output.PosRaw} > {output.SnpPos}
        """

rule LeflerSnpsMatrixEqtl:
    input:
        GenotypeRaw = "GTEX_renalysis/LeflerSnps/Genotypes.012",
        IndvRaw = "GTEX_renalysis/LeflerSnps/Genotypes.012.indv",
        PosRaw = "GTEX_renalysis/LeflerSnps/Genotypes.012.pos",
        GenotypeFinished = "GTEX_renalysis/LeflerSnps/MatrixEQTL/Genotypes.txt",
        SnpPos = "GTEX_renalysis/LeflerSnps/MatrixEQTL/snploc.txt",
        expression = "GTEX_renalysis/LeflerSnps/MatrixEQTL/expression.txt",
        geneloc = "GTEX_renalysis/LeflerSnps/MatrixEQTL/geneloc.txt",
        GenotypeRaw_m = "GTEX_renalysis/LeflerSnps/Genotypes.Control.012",
        IndvRaw_m = "GTEX_renalysis/LeflerSnps/Genotypes.Control.012.indv",
        PosRaw_m = "GTEX_renalysis/LeflerSnps/Genotypes.Control.012.pos",
        GenotypeFinished_m = "GTEX_renalysis/LeflerSnps/MatrixEQTL/Genotypes.control.txt",
        SnpPos_m = "GTEX_renalysis/LeflerSnps/MatrixEQTL/snploc.control.txt",
        SNP_file_name = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedSnps.PanTro5.snps",
        snps_location_file_name = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedSnps.PanTro5.snploc",
        expression_file_name = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt",
        gene_location_file_name = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.geneloc.txt",
        covariates_file_name = "../../output/Covariates/0GenotypePCs_and_10RNASeqPCs.covariates",
        errorCovariance_file = "eQTL_mapping/Kinship/GRM.cXX.txt",
        SNP_file_name_matched = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedMatchedSnps.PanTro5.snps",
        snps_location_file_name_matched = "eQTL_mapping/SharedPolymorphisms/SpeciesSharedMatchedSnps.PanTro5.snploc"
    output:
        "../../output/LeflerTestedSnps.Chimp.cis.tsv"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/MatrixEqtl_CisAndTrans.AllPvals.Chimp.R
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/MatrixEqtl_CisAndTrans.AllPvals.Human.R
        """

rule CopySnpPosKeyHumanLeflerSnps:
    input:
        "GTEX_renalysis/LeflerSnps/MatrixEQTL/snploc.txt",
    output:
        "../../output/LeflerTestedSnps.Human.snploc.txt",
    shell:
        "cp {input} {output}"

rule GetGenotypesForSharedSnpBoxplot:
    input:
        GtexVcf = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz",
        GtexVcfTbi = "GTEX_renalysis/vcf/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz.tbi",
        covariates = "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Heart_Left_Ventricle.v8.covariates.txt",
        ChimpVcf = "PopulationSubstructure/ReferencePanelMerged.annotated.vcf.gz",
        expression = "GTEX_renalysis/LeflerSnps/MatrixEQTL/expression.txt",
    params:
        human_snp = "chr5:128978923",
        chimp_snp = "5:128875253",
    output:
         human_gt = "../../output/SharedSnpExampleBoxplot.human.012",
         chimp_gt = "../../output/SharedSnpExampleBoxplot.chimp.012",
         expression = "../../output/SharedSnpExampleBoxplot.human.expression.txt"
    shell:
        """
        bcftools view -S <(cat {input.covariates} | head -1 |transpose | awk 'NR>1')  -r {params.human_snp} {input.GtexVcf} | vcftools --vcf - --012 --out ../../output/SharedSnpExampleBoxplot.human
        bcftools view -S <(bcftools query -l {input.ChimpVcf} | grep "ThisStudy" | grep -v "MD_And") -r {params.chimp_snp} {input.ChimpVcf} | vcftools --vcf - --012 --out ../../output/SharedSnpExampleBoxplot.chimp
        cat {input.expression} | awk 'NR==1 || $1=="ENSG00000113396.12"' | transpose > {output.expression}
        """

# rule DownloadGTExSummaryStatsAllTissues:
#     output:
#         expand()
#     log:
#         "logs/DownloadGTExSummaryStats.log"
#     shell:
#         """
#         wget ...*
#         """


# rule MergeGTExSummaryStatsAllTissues:
#     input:
#         expand()
#     log:
#         "logs/CombineGTExSummaryStats.log"
#     output:
#         ""
#     shell:
#         """
#         """
