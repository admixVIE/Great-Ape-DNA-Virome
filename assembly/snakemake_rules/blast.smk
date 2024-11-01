rule blast_fragments:
    input:
        contig="../results/{experiment}/05_genes/{sample}.genes.fasta"
    output:
        hits="../results/{experiment}/06_blast_hits/{sample}.hits",
        names="../results/{experiment}/06_blast_hits/{sample}.names"
    log:
        "log/blast_{experiment}_{sample}.log"
    threads: 5
    resources: mem_mb=8000
    params:
        blastdb="../db/blastdb/ref_viruses_rep_genomes"
    shell:
        """
        module load ncbiblastplus edirect
        #Check if file is empty
        if [ -s {input.contig} ]; then

        blastn -query {input.contig} -db {params.blastdb} -num_threads {threads} -outfmt 7 -out {output.hits}

        #Extract best hits
        head -4 {output.hits} | tail -1 | cut -f 2 -d ":" | tr -d " " | tr "," "\t" | awk '{{$1 = "contig_id"; $2 = "acc_number"; print}}' | tr " " "\t" > $TMPDIR/besthits.tab
        grep -A 5 "^# BLASTN" {output.hits} | grep -v "^#" | grep -v "^-" >> $TMPDIR/besthits.tab || true
        tail -n +2 $TMPDIR/besthits.tab | cut -f 2 > $TMPDIR/accs_tmp
        sort -u $TMPDIR/accs_tmp > $TMPDIR/accs

        #Get the names of the best hits
        echo -e "acc_number\t{wildcards.sample}" > {output.names}

        
        #skip if there is nothing in accs
        if [ -s $TMPDIR/accs ]; then
        for i in $(cat $TMPDIR/accs)
        do
        echo -e "$i\t$(echo $i | epost -db nuccore | efetch -format xml | xtract -pattern Bioseq-set_descr -element Org-ref_taxname | awk '{{if ($i =="") print "unknown"; else print}}')" >> {output.names}
        done
    else
        touch {output.hits} {output.names}
        fi

    else
        touch {output.hits} {output.names}
        fi

        """
