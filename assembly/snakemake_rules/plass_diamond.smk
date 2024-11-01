rule plass_diamond:
    input:
        assembly="../results/{experiment}/04b_plass/{sample}.faa",
        lca="../results/{experiment}/05b_plass_taxonomy/{sample}_lca.tab"
    output:
        filtered="../results/{experiment}/05c_plass_filtered/{sample}.faa",
        dmnd="../results/{experiment}/05c_plass_filtered/{sample}.dmnd"
    threads: 8
    resources:
        mem_mb=16000,
        time="2-00:00:00"
    params:
        filter="bin/filter_taxonomy.py",
        filter_matches="bin/filter_matches.py",
        nr="/lisc//scratch/cube/Viruses_of_great_apes/db/nr_diamond/nr.dmnd",
        acc2tax="/lisc/scratch/cube/Viruses_of_great_apes/db/acc2prot/prot.accession2taxid.edit2"
    shell:
        """

        #Check if the input file is empty
        if [ ! -s {input.lca} ]; then
            touch {output.filtered}
            exit 0
        fi

        module unload python3
        module load ete3

        #Filter out the undesired proteins
        python3 {params.filter} --infile {input.lca} --taxcol 2 --taxidsL 10239 --taxidsH 1 --unclassified yes | cut -f 1 > $TMPDIR/good_names

        #Take the subset of selected genes
        module load seqtk
        seqtk subseq {input.assembly} $TMPDIR/good_names > $TMPDIR/good_contigs.faa

        #Diamond filter the subset of contigs
        #module load conda
        #conda activate diamond
        module load diamond
        diamond blastp -d {params.nr}\
                       -q $TMPDIR/good_contigs.faa\
                       --outfmt 6 qseqid qlen sseqid slen evalue bitscore length staxids \
                       --max-target-seqs 1\
                       -o $TMPDIR/matches1.tsv\
                       --threads {threads}\
                       --masking 0
        cp $TMPDIR/matches1.tsv {output.dmnd}
        #if [ ! -s $TMPDIR/matches1.tsv ]; then
        #    cat $TMPDIR/matches1.tsv > {output.dmnd}
        #else
        #    touch {output.dmnd}
        #fi

       #Exit if there are no results
       if [ ! -s $TMPDIR/matches1.tsv ]; then
           touch {output.filtered}
           exit 0
       fi

        python3 {params.filter_matches} $TMPDIR/matches1.tsv | awk '{{OFS="\t"; split($3,a,".");$3=a[1]".";print}}' | cut -f 1-7 | sort -k3b,3 > $TMPDIR/dmnd.filter.1

        #Exit if there are no results
        if [ ! -s $TMPDIR/dmnd.filter.1 ]; then
            touch {output.filtered}
            exit 0
        fi


        cut -f 3 $TMPDIR/dmnd.filter.1 > $TMPDIR/accessions
        grep -F -f $TMPDIR/accessions {params.acc2tax} > $TMPDIR/taxids.txt
        sort -k1b,1 $TMPDIR/taxids.txt > $TMPDIR/taxids
        join -1 3 -2 1 -a1 $TMPDIR/dmnd.filter.1 $TMPDIR/taxids > $TMPDIR/dmnd.filter.2

        python3 {params.filter_matches} $TMPDIR/dmnd.filter.2 > $TMPDIR/good_contigs.dmnd.filter.tsv

        cut -f 2 -d " " $TMPDIR/good_contigs.dmnd.filter.tsv > $TMPDIR/names_diamond

        grep "^>" $TMPDIR/good_contigs.faa | cut -f 2 -d ">" | cut -f 1 -d " " | sort > $TMPDIR/all_proteins

        cut -f 1 $TMPDIR/matches1.tsv | sort > $TMPDIR/all_diamond

        comm -23 $TMPDIR/all_proteins $TMPDIR/all_diamond > $TMPDIR/no_diamond

        cat $TMPDIR/names_diamond $TMPDIR/no_diamond > $TMPDIR/names
        seqtk subseq {input.assembly} $TMPDIR/names > {output.filtered}

        rm -r $TMPDIR
        """
