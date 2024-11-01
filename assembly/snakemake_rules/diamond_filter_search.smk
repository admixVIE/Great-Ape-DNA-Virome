rule diamond_search:
    input:
        proteins="../results/{experiment}/05_genes/{sample}.proteins.fasta",
        genes="../results/{experiment}/05_genes/{sample}.genes.fasta"
    output:
        diamondres="../results/{experiment}/05a_filtered_genes/{sample}.dmnd.tsv"
        #temp(diamondres="../results/{experiment}/05a_filtered_genes/{sample}.dmnd.tsv")
    log:
        "log/diamond_{experiment}_{sample}.log"
    threads: 4
    resources: mem_mb=16000,
               time="2-00:00:00"
    params:
        filtermatches="bin/filter_matches.py",
        ac2tax="../db/acc2prot/prot.accession2taxid.edit",
        diamonddb="../db/nr_diamond/nr.dmnd"
    shell:
        """
        module unload python3
        module load seqtk ete3

        #Check if the input file is empty
        if [ ! -s {input.proteins} ]; then
            touch {output.diamondres}
            exit 0
        fi

        #Diamond search agains NCBI nr database
        #Copy db to tmpdir
        cp {params.diamonddb} $TMPDIR/db.dmnd

        module load diamond
        #Protein search
        diamond blastp -d $TMPDIR/db.dmnd \
                          -q {input.proteins} \
                          --outfmt 6 qseqid qlen sseqid slen evalue bitscore length staxids \
                          --max-target-seqs 1 \
                          -o $TMPDIR/matches1.tsv \
                          --masking 0



        #Extract genes that did not have a match using blastp
        cut -f1 $TMPDIR/matches1.tsv | sort > $TMPDIR/dmnd.yes
        grep "^>" {input.genes} | cut -f1 -d " "| tr -d ">" | sort > $TMPDIR/all.genes
        comm -13 $TMPDIR/dmnd.yes $TMPDIR/all.genes > $TMPDIR/dmnd.no
        seqtk subseq {input.genes} $TMPDIR/dmnd.no > $TMPDIR/dmnd_no.fasta



        #Run diamond blastx on the genes not having hits with blastp
        diamond blastx -d $TMPDIR/db.dmnd \
                       -q $TMPDIR/dmnd_no.fasta \
                       --outfmt 6 qseqid qlen sseqid slen evalue bitscore length staxids \
                       --max-target-seqs 1 \
                       -o $TMPDIR/matches2.tsv \
                       --masking 0

        


        #Merge results from two diamond searches
        cat $TMPDIR/matches1.tsv $TMPDIR/matches2.tsv | sort > $TMPDIR/matches.tsv

        #log information
        echo "the number of lines in a file with blastx results"
        echo $(wc -l $TMPDIR/matches2.tsv)

        #Copy the result to an appropriate folder
        mv $TMPDIR/matches.tsv {output.diamondres}

        """
