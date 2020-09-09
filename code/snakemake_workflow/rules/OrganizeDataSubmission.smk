rule CopyNovelRNASeqData:
    input:
        lambda wildcards: GEORNASeqBasenameToFastq[wildcards.basename]
    output:
        "DataForSubmission/GEO/Fastq/{basename}"
    shell:
        """
        cp {input} {output}
        """

rule CheckSumRNASeq:
    input:
        RNASeqFastqList = expand( "DataForSubmission/GEO/Fastq/{basename}", basename=GEORNASeqBasenameToFastq.keys())
    output:
        "DataForSubmission/GEO/FastqChecksums.txt"
    shell:
        """
        md5sum {input} > {output}
        """

rule MetaDataFormatHelper:
    input:
        "DataForSubmission/GEO/FastqChecksums.txt",
        "../../data/NovelRNASeqFiles.tsv"
    output:
        "DataForSubmission/GEO/MetadataHelper.tsv"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript 
        """

rule CopyProcessedFiles:
    input:
        ChimpTable = "PowerAnalysis/Subread/Chimp.subread.txt.gz",
        HumanTable = "PowerAnalysis/Subread/Human.subread.txt.gz",
        ChimpRawTable = "RNASeq/STAR/CountTable.txt.gz",
        ChimpNormalizedTable = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.phenotypes.txt"
    output:
        ChimpTable = "DataForSubmission/GEO/OrthologousExonRawGeneCounts_Chimpanzee.tsv.gz",
        HumanTable = "DataForSubmission/GEO/OrthologousExonRawGeneCounts_Human.tsv.gz",
        ChimpRawTable = "DataForSubmission/GEO/ChimpanzeesRawGeneCount_eQTL_Mapping.tsv.gz",
        ChimpNormalizedTable = "DataForSubmission/GEO/ChimpanzeesQQNorm_eQTL_Mapping.tsv.gz"
    shell:
        """
        cp {input.ChimpTable} {output.ChimpTable}
        cp {input.HumanTable} {output.HumanTable}
        cp {input.ChimpRawTable} {output.ChimpRawTable}
        cat {input.ChimpNormalizedTable} | gzip - > {output.ChimpNormalizedTable}
        """

rule GetEVA_Approved_Reference:
    """
    I called variants using a reference with non-chromosomal contigs deleted.
    This is probably not best practice, and not approved by European Variation
    Archive for submission. I will need to update the vcf header to include all
    contigs from the ref genome GCA_000001515.5
    """
    output:
        fa = "Misc/EnsemblRef/Ref.fa",
        faidx =  "Misc/EnsemblRef/Ref.fa.fai",
        ContigInfo = "Misc/EnsemblRef/Ref.report.tab"
    shell:
        """
        wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/515/GCA_000001515.5_Pan_tro_3.0/GCA_000001515.5_Pan_tro_3.0_genomic.fna.gz | gunzip > {output.fa}
        samtools faidx {output.fa}
        wget -O {output.ContigInfo} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/515/GCA_000001515.5_Pan_tro_3.0/GCA_000001515.5_Pan_tro_3.0_assembly_report.txt
        """

rule CreateChromConversionTable:
    input:
        fai =config["ref"]["genome"] + ".fai",
        ContigInfo = "Misc/EnsemblRef/Ref.report.tab"
    output:
        "MiscOutput/ChrConversionKey.tab"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/CheckRefSeqs.R {input.fai}
        """

rule RenameChromsInVcf:
    input:
        filtered= "PopulationSubstructure/ReferencePanelMerged.annotated.vcf.gz",
        faidx =  "Misc/EnsemblRef/Ref.fa.fai",
        ConversionKey = "MiscOutput/ChrConversionKey.tab"
    output:
        vcf = "DataForSubmission/EVA/ToSubmit.vcf.gz",
        tbi  = "DataForSubmission/EVA/ToSubmit.vcf.gz.tbi",
        NewHeader =  "DataForSubmission/EVA/ToSubmit.Header.vcf"
    shell:
        """
        cat <(bcftools view -h {input.filtered} | head -5) <(awk -F'\\t' '{{ print "##contig=<ID="$1",length="$2">"  }}' {input.faidx}) <(bcftools view -h {input.filtered} | sed -n '/##phasing=/,$p') > {output.NewHeader}
        bcftools annotate  --rename-chrs {input.ConversionKey} {input.filtered} | bcftools reheader -h {output.NewHeader} - | bcftools view -O z  -S <(bcftools query -l {input.filtered} | grep 'ThisStudy') > {output.vcf}
        tabix -p vcf {output.vcf}
        """


