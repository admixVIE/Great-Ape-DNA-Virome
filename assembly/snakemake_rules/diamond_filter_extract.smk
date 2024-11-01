rule extract_viral_proteins:
    input:
        proteins="../results/{experiment}/05_genes/{sample}.proteins.fasta",
        diamond="../results/{experiment}/05a_filtered_genes/{sample}.dmnd.filter.tsv",
        diamondall="../results/{experiment}/05a_filtered_genes/{sample}.dmnd.tsv"
    output:
        outprot="../results/{experiment}/05a_filtered_genes/{sample}.proteins.fasta"
    log:
        "log/diamond_filter_tmp_{experiment}_{sample}.log"
    threads: 1
    resources: mem_mb=2000,
               time="5:00:00"
    params:
        filtermatches="bin/filter_matches.py",
        ac2tax="../db/acc2prot/prot.accession2taxid.edit2"
    shell:
        """
        #Check if the input file is empty
        if [ ! -s {input.diamond} ]; then
            touch {output.outprot}
            exit 0
        fi

        module load seqtk

        #Hits that are not non-viral and have diamond hits
        cut -f 2 -d " " {input.diamond} > $TMPDIR/names_diamond

        #all proteins
        grep "^>" {input.proteins} | cut -f 2 -d ">" | cut -f 1 -d " " | sort > $TMPDIR/all_proteins

        #Proteins that have diamond hits
        cut -f 1 {input.diamondall} | sort >$TMPDIR/all_diamond

        #Proteins that don't have diamond hits
        comm -23 $TMPDIR/all_proteins $TMPDIR/all_diamond > $TMPDIR/no_diamond

        #Combined proteins that don't have diamond hits and those that were not eliminated by diamond
        cat $TMPDIR/names_diamond $TMPDIR/no_diamond > $TMPDIR/names

        #Extract the desired proteins in fasta format
        seqtk subseq {input.proteins} $TMPDIR/names > {output.outprot}
        """
