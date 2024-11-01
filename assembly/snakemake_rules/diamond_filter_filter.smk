rule diamond_filter:
    input:
        proteins="../results/{experiment}/05_genes/{sample}.proteins.fasta",
        diamond="../results/{experiment}/05a_filtered_genes/{sample}.dmnd.tsv"
    output:
        diamondres="../results/{experiment}/05a_filtered_genes/{sample}.dmnd.filter.tsv"
    log:
        "log/diamond_filter_{experiment}_{sample}.log"
    threads: 1
    resources: mem_mb=4000,
               time="05:00:00"
    params:
        filtermatches="bin/filter_matches.py",
        ac2tax="../db/acc2prot/prot.accession2taxid.edit2"
    shell:
        """
        #Check if the input file is empty
        if [ ! -s {input.diamond} ]; then
            touch {output.diamondres}
            exit 0
        fi

        #Remove all proteins that are probably not viral, round 1
        module unload python3
        module load ete3

        #Filter with python script, remove the taxid column
        python3 {params.filtermatches} {input.diamond} | awk '{{OFS="\t"; split($3,a,"."); $3=a[1]".";print}}' | cut -f 1-7 | sort -k3,3 > $TMPDIR/dmnd.filter.1

        #Extract accession number without version from what is left after first filtering run
        cut -f 3 $TMPDIR/dmnd.filter.1  > $TMPDIR/accessions

        #Grep search for taxids
        grep -F -f $TMPDIR/accessions {params.ac2tax} > $TMPDIR/taxids.txt || true
        touch $TMPDIR/taxids.txt
        sort -k1 $TMPDIR/taxids.txt > $TMPDIR/taxids

        #Merge diamond results with taxids
        join -1 3 -2 1 -a1 $TMPDIR/dmnd.filter.1 $TMPDIR/taxids > $TMPDIR/dmnd.filter.2

        #Remove all proteins that are probably not viral, round 2
        python3 {params.filtermatches} $TMPDIR/dmnd.filter.2 > {output.diamondres}
        """
