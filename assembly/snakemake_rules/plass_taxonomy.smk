rule plass_taxonomy:
    input:
        assembly="../results/{experiment}/04b_plass/{sample}.faa"
    output:
        lca="../results/{experiment}/05b_plass_taxonomy/{sample}_lca.tab"
    threads:16
    resources:
        mem_mb=400000,
        time="4-20:00:00"
    params:
        nr="/lisc/scratch/cube/Viruses_of_great_apes/db/nr_mmseqs/nr"
    shell:
        """
        module load mmseqs2

        mkdir $TMPDIR/res $TMPDIR/tmp $TMPDIR/inp

        module load conda
        conda activate seqkit

        seqkit rmdup -s < {input.assembly} > $TMPDIR/inp/inp.faa

        mmseqs easy-taxonomy $TMPDIR/inp/inp.faa {params.nr} $TMPDIR/res/res $TMPDIR/tmp --threads {threads} --split-memory-limit 400G

        cp $TMPDIR/res/res_lca.tsv {output.lca}

        rm -r $TMPDIR
        """
