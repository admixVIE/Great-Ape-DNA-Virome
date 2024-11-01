rule diamond_coassembly:
    input:
        contigs="../results/{experiment}/04_contigs/coassembly.fasta",
        proteins="../results/{experiment}/05_genes/coassembly.proteins.fasta"
    output:
        dmnd=temp("../results/{experiment}/05_coassembly/coassembly.dmnd")
        #annotated_contigs="../results/{experiment}/05_coassembly/annotated_contigs.tab"
    log:
        "log/coassembly_dmnd_{experiment}.log"
    threads: 8
    params:
        diamond_db="../db/nr_diamond/nr.dmnd",
        ac2tax="../db/acc2prot/prot.accession2taxid.edit2",
        fasta_to_dmnd="../bin/fasta_to_dmndlike.py",
        annotate="../bin/annotate_contigs.py"
    resources:
        mem_mb=50000,
        time="1-00:00:00"
    shell:
        """
        module load conda

        conda activate diamond-2.1.6

        diamond blastp -d {params.diamond_db} \
                -q {input.proteins} \
                --outfmt 6 qseqid qlen sseqid slen evalue bitscore length staxids \
                --max-target-seqs 1 \
                -o {output.dmnd} \
                --masking 0 \
                -p {threads}
        """

rule diamond_contigs:
    input:
        proteins="../results/{experiment}/05_genes/coassembly.proteins.fasta",
        dmnd="../results/{experiment}/05_coassembly/coassembly.dmnd"
    output:
        annotated_contigs="../results/{experiment}/05_coassembly/annotated_contigs.tab"
    log:
        "log/coassembly_contigs_{experiment}.log"
    threads: 1
    params:
        ac2tax="../db/acc2prot/prot.accession2taxid.edit2",
        fasta_to_dmnd="bin/fasta_to_dmndlike.py",
        annotate="bin/annotate_contigs.py"
    resources:
        mem_mb=5000,
        time="1-00:00:00"
    shell:
        """
        ##Improve matches

        #Remove the version number of the genbank entry
        cat {input.dmnd} | awk '{{OFS="\t"; split($3,a,"."); $3=a[1]"."; print}}' | cut -f 1-7 | sort -k3,3 > $TMPDIR/matches_edit1.tsv

        #Extract the accession without the version
        cut -f 3 $TMPDIR/matches_edit1.tsv > $TMPDIR/accessions

        #Search for taxids
        grep -F -f $TMPDIR/accessions {params.ac2tax} > $TMPDIR/taxids.txt || true
        sort -k 1b,1 $TMPDIR/taxids.txt > $TMPDIR/taxids

        #Merge diamond matches with the extended list of taxids
        join -1 3 -2 1 -a1 $TMPDIR/matches_edit1.tsv $TMPDIR/taxids | awk 'OFS="\t" {{j=$1; d=$2; t=$3; $1=d; $2=t; $3=j;print}}' > $TMPDIR/matches_extended.tsv

        #########################################################################################################################################################

        module load ete3
        cut -f1 $TMPDIR/matches_extended.tsv | sort > $TMPDIR/dmnd.yes

        grep "^>" {input.proteins} | cut -f 1 -d " " | tr -d ">" | sort > $TMPDIR/all.prots

        comm -13 $TMPDIR/dmnd.yes $TMPDIR/all.prots > $TMPDIR/dmnd.no

        python3 {params.fasta_to_dmnd} $TMPDIR/dmnd.no > $TMPDIR/dmnd.no.dmnd

        cat $TMPDIR/dmnd.no.dmnd $TMPDIR/matches_extended.tsv | sort > $TMPDIR/all.dmnd

        python3 {params.annotate} $TMPDIR/all.dmnd > {output.annotated_contigs}
        """
