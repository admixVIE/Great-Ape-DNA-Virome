rule plass_mapVOG:
    input:
        filtered="../results/{experiment}/05c_plass_filtered/{sample}.faa"
    output:
        vog="../results/{experiment}/06b_plass_vog/{sample}.m8",
        best="../results/{experiment}/06b_plass_vog/{sample}_best.tsv"
    threads: 16
    resources:
        mem_mb=400000,
        time="10:00:00"
    params:
        vog_mmseqs="/lisc/scratch/cube/Viruses_of_great_apes/db/vogdb_mmseqs/vogdb",
        lca="/lisc/scratch/cube/Viruses_of_great_apes/db/vogdb_mmseqs/vog.lca.tsv",
        extract_best="/lisc/scratch/cube/Viruses_of_great_apes/scripts/bin/extract_best_m8.py",
        merge="/lisc/scratch/cube/Viruses_of_great_apes/scripts/bin/join_lca.R"
    shell:
        """
        #Check if the input file is empty
        if [ ! -s {input.filtered} ]; then
            touch {output.vog} {output.best}
            exit 0
        fi

        module load mmseqs2

        mkdir $TMPDIR/tmp

        #Map proteins to vogs
        mmseqs easy-search {input.filtered} {params.vog_mmseqs} {output.vog} $TMPDIR/tmp --threads {threads} --split-memory-limit 200G

        if [ ! -s {output.vog} ]; then
            touch {output.vog} {output.best}
            exit 0
        fi

        #Process the resulting file
        python3 {params.extract_best} --infile {output.vog} --lca {params.lca} > $TMPDIR/best_hits

        module load R
        Rscript {params.merge} {params.lca} $TMPDIR/best_hits {output.best} #Instead of best_hits it was output_vog

        rm -r $TMPDIR
        """
