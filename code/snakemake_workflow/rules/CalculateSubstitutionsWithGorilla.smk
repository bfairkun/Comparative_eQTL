rule GetAppris:
    output:
        "../../data/Appris.principal.isoforms.txt"
    shell:
        """
        wget -O {output} "http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt"
        """

rule AlignCDS_andCountSubstitutions:
    """
    Rscript to download coding sequences from Biomart and align them with DECIPHER
    library and use a system command from within the Rscript to the SNAP.pl
    script to count the subtitutions and write to a summary file for each gene.
    Then concatentate the summary files into a single output.
    """
    input:
        appris = "../../data/Appris.principal.isoforms.txt"
    output:
        "../../output/GorillaChimpHumanFixedSubstitutionsCount.gz"
    log:
        "logs/AlignCDS_andCountSubstitutions.log"
    shell:
        """
        mkdir -p CDS_alignments/ CDS_alignment_SNAPfiles/
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/AlignCDS.ChimpHumanGorilla.R &> {log}
        find ./CDS_alignment_SNAPfiles/summary.* | xargs awk 'FNR>1 && FNR<5 {{ print $0, FILENAME }}' | gzip - > {output}
        """
