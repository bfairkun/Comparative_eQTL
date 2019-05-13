rule STAR_to_leafcutter_junc:
    input:
        "RNASeq/STAR/{sample}/SJ.out.tab"
    output:
        "sQTL_mapping/juncfiles/{sample}.junc"
    log:
        "logs/sQTL_mapping/STAR_to_leafcutter_junc/{sample}.log"
    shell:
        """
        my awk code
        """

rule make_leafcutter_juncfile:
    input:
        expand("")

rule leafcutter_cluster:

rule leafcutter_prepare_phenotype_table:

rule preapre_MatrixEQTL_for_sQTL:

rule MatrixEQTL_for_sQTL:
