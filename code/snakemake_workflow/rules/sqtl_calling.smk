# include: "common.smk"
# rule all:
#     input:
#         expand ("sQTL_mapping/juncfiles/{sample}.junc", sample=samples["sample"] ),
#         "sQTL_mapping/leafcutter/juncfilelist.txt",
#         "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz",

rule STAR_to_leafcutter_junc:
    input:
        "RNASeq/STAR/{sample}/SJ.out.tab"
    output:
        "sQTL_mapping/juncfiles/{sample}.junc"
    log:
        "logs/sQTL_mapping/STAR_to_leafcutter_junc/{sample}.log"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$4==1 && $1!="MT" {{ print $1,$2,$3,".",$7,"+" }} $4==2&& $1!="MT" {{ print $1,$2,$3,".",$7,"-" }}' {input} > {output}
        """

rule make_leafcutter_juncfile:
    input:
        expand ("sQTL_mapping/juncfiles/{sample}.junc", sample=samples["sample"] ),
    output:
        "sQTL_mapping/leafcutter/juncfilelist.txt"
    params:
        SamplesToRemove = "MD_And"
    run:
        import os
        with open(output[0], "w") as out: 
            for filepath in input:
                samplename = os.path.basename(filepath).split(".junc")[0]
                if samplename != params.SamplesToRemove:
                    out.write(filepath + '\n')

rule leafcutter_cluster:
    input:
        "sQTL_mapping/leafcutter/juncfilelist.txt",
    output:
        "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz",
        "sQTL_mapping/leafcutter/clustering/leafcutter_perind_numers.counts.gz"
    log:
        "logs/sQTL_mapping/leafcutter_cluster.log"
    shell:
        """
        leafcutter_cluster.py -j {input} -r sQTL_mapping/leafcutter/clustering/ &> {log}
        """

rule leafcutter_prepare_phenotype_table:
    input:
        counts = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz",
        blacklist_chromosomes = "scratch/chromsomseblacklist.txt"
    output:
        "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_Catted.txt"
    shell:
        """
        ~/miniconda3/bin/python2.7 ~/CurrentProjects/leafcutter/scripts/prepare_phenotype_table.py -p 13 --ChromosomeBlackList {input.blacklist_chromosomes}  {input.counts}

        cat <(head -1 sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_chr3) <(awk 'FNR>1' sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_chr*) > 
        """

rule preapre_MatrixEQTL_for_sQTL:
    input:
        "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_Catted.txt"
    output:
        phenotypes = "sQTL_mapping/MatrixEQTL/sQTL_phenoTable.txt",
        intron_locs = "sQTL_mapping/MatrixEQTL/sQTL_intron.locs"


# rule MatrixEQTL_for_sQTL:
