rule GetAppris:
    output:
        "../../data/Appris.principal.isoforms.txt"
    shell:
        """
        wget -O {output} "http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt"
        """

rule AlignCDS:
    """
    Rscript to download coding sequences from Biomart. Align them with DECIPHER library
    """
    input:
        appris = "../../data/Appris.principal.isoforms.txt"
