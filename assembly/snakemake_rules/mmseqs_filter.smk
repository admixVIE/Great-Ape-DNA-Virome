rule mmseqs_filter:
    input:
        contigs="../results/{experiment}/04_contigs/{sample}.fasta",
    output:
        lca="../results/{experiment}/045_mmseq_taxonomy/{sample}_lca.tsv",
        good_contigs="../results/{experiment}/046_filtered_contigs/{sample}.fasta"
    threads:8
    resources:
        mem_mb=230000,
        time="2-10:00:00"
    params:
        refdb="/lisc/scratch/cube/Viruses_of_great_apes/db/nr_mmseqs/nr",
        extractscript="bin/filter_taxonomy.py"
    shell:
        """
        module load mmseqs2

        mkdir $TMPDIR/tmp

        #Taxonomy determination using mmseqs
        /usr/bin/time mmseqs easy-taxonomy {input.contigs} {params.refdb} $TMPDIR/{wildcards.sample} $TMPDIR/tmp --threads {threads} --split-memory-limit 230G || touch {output.lca} {output.good_contigs} #; exit 0

        cp $TMPDIR/{wildcards.sample}* $(dirname {output.lca})

        #Extract only the good contigs
        #awk -F'\t' '$4 == "root"' {output.lca} | cut -f 1 > $TMPDIR/good_names
        #awk -F'\t' '$4 == "unclassified"' {output.lca} | cut -f 1 >> $TMPDIR/good_names
        module unload python3
        module load ete3
        python3 {params.extractscript} --infile {output.lca} --taxcol 2 --taxidsL 10239 --taxidsH 1 --unclassified yes | cut -f 1 > $TMPDIR/good_names

        module load seqtk
        seqtk subseq {input.contigs} $TMPDIR/good_names > $TMPDIR/contigs.fasta
        cp $TMPDIR/contigs.fasta {output.good_contigs}

        ls -lah $TMPDIR
        #rm -r $TMPDIR
        """
