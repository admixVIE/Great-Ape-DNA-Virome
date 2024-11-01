rule filter_unmapped:
    input:
        fasta="../results/{experiment}/04_unmapped/{sample}.cdhit.fasta.gz"
    output:
        taxonomy="../results/04_unmapped_taxonomy/{sample}.lca",
        filtered="../results/04_unmapped_filtered/{sample}.fasta.gz"
    threads:8
    resources:
        mem_mb=120000,
        time="04_00:00:00"
    params:
        refdb="/lisc/scratch/cube/Viruses_of_great_apes/db/nr_mmseqs/nr",
        extractscript="bin/filter_taxonomy.py"
    shell:
        """
        module load mmseqs2

        mkdir $TMPDIR/tmp

        #search the reference database
        mmseqs easy-taxonomy {input.fasta} {params.refdb} $TMPDIR/{wildcards.sample} $TMPDIR/tmp --threads {threads} --split-memory-limit 120G

        cp $TMPDIR/{wildcards.sample}_lca.tsv {output.taxonomy}

        #Extract good reads
        python3 {params.extractscript} --infile $TMPDIR/{wildcards.sample}_lca.tsv --taxcol 2 --taxidsL 10239 --taxidsH 1 --unclassified yes | cut -f 1 > $TMPDIR/good_names

        module load seqtk
        seqtk subseq {input.fasta} $TMPDIR/good_names > $TMPDIR/good.fasta
        gzip $TMPDIR/good.fasta

        cp $TMPDIR/good.fasta.gz {output.filtered}

        rm -r $TMPDIR

        """
