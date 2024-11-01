rule plass_desired_taxa:
    input:
        vog="../results/{experiment}/06b_plass_vog/{sample}_best.tsv",
        mmseqs="../results/{experiment}/05b_plass_taxonomy/{sample}_lca.tab"
    output:
        vogs="../results/{experiment}/06c_selected/{sample}_vog.tsv",
        mmseqs="../results/{experiment}/06c_selected/{sample}_mmseqs.tsv"
    threads:1
    resources:
        mem_mb=2000,
        time="02:00:00"
    params:
        select="bin/desired_taxonomies.py",
        expand="bin/expand_lineage.py"
    shell:
        """
        module unload python3
        module load ete3

        #Taxids
        #Adenoviridae 10508
        #Circoviridae 39724
        #Hepadnaviridae 10404
        #Orthoherpesviridae 3044472
        #Papillomaviridae 151340
        #Parvoviridae 10780
        #Poxviridae 10240
        #Retroviridae 11632

        #Exit if there are no results
        if [ ! -s {input.vog} ]; then
            touch {output.vogs} {output.mmseqs}
            exit 0
        fi

        #Exit if there are no results
        if [ ! -s {input.mmseqs} ]; then
            touch {output.vogs} {output.mmseqs}
            exit 0
        fi

        python3 {params.select} --infile {input.vog} --lcacol 5 --taxids 10508,39724,10404,3044472,151340,10780,10240,11632 --head yes > {output.vogs}

        python3 {params.select} --infile {input.mmseqs} --lcacol 2 --taxids 10508,39724,10404,3044472,151340,10780,10240,11632 > $TMPDIR/tmp

        python3 {params.expand} --infile $TMPDIR/tmp --taxcol 2 --head no > {output.mmseqs}

        rm -r $TMPDIR
        """
