rule contig_filtering:
    input:
        proteins="../results/{experiment}/05_genes/{sample}.proteins.fasta",
        diamond="../results/{experiment}/05_genes/{sample}.proteins.fasta"
    output:
        canno="../results/{experiment}/05b_contig_filter/{sample}.contig.anno"
    log:
        "log/contig_anno_{sample}.log"
    threads: 1
    resources: mem_mb=2000
    params:
        anno_contigs="bin/annotate_contigs.py",
        fastadmnd="bin/fasta_to_dmndlike.py"
    shell:
        """
        #Check if the input file is empty
        if [ ! -s {input.diamond} ]; then
            touch {output.canno}
            exit 0
        fi

        #Diamond proteins
        cut -f1 {input.diamond} | sort > $TMPDIR/dmnd.yes

        #All proteins
        grep "^>" {input.proteins} | cut -f 1 -d " "| tr -d ">" | sort > $TMPDIR/all.prots

        #Proteins that don't have diamond hits
        comm -13 $TMPDIR/dmnd.yes $TMPDIR/all.prots > $TMPDIR/dmnd.no

        #Prepare format
        python3 {params.fastadmnd} $TMPDIR/dmnd.no > $TMPDIR/dmnd.no.dmnd

        #Merge all dmnd outputs
        cat $TMPDIR/dmnd.no.dmnd {input.diamond} | sort > $TMPDIR/all.dmnd

        #Annotate contigs
        python3 {params.anno_contigs} $TMPDIR/all.dmnd > {output.canno}
        """
