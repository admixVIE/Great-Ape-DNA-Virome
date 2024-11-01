rule prepare:
    input:
        viruses="../results/{experiment}/06_vog_map/{sample}.voglin",
        abund="../results/{experiment}/05a_gene_abundances/{sample}.abundances.tab",
        voglist="../results/{experiment}/06_vog_map/{sample}.voglist"
    output:
        table="../results/{experiment}/07_tables/viruses/{sample}.tab"
    log:
        "log/prepare_{experiment}_{sample}.log"
    threads: 1
    resources: mem_mb=1000,
               time="3:00:00"
    params:
        arrange="bin/table_from_virlin.py",
        phageterms="etc/phage_lineage_terms.txt"
    shell:
        """
        #Check fi the file is not empty
        if ! [ -s {input.viruses} ];then
            touch {output.table}
            exit 0
        fi

        #Filter out phages
        grep -v -f {params.phageterms} {input.viruses} | sort -u > $TMPDIR/viruses || true

        #Prepare and join voglist and abundance files
        sort -k1,1 {input.abund} > $TMPDIR/abund.sorted
        sort -k2,2 {input.voglist} > $TMPDIR/voglist.sorted
        join -1 1 -2 2 $TMPDIR/abund.sorted $TMPDIR/voglist.sorted -t $'\t' > $TMPDIR/contig_abund
        sort -k1,1 $TMPDIR/contig_abund > $TMPDIR/contig_abund.sorted

        #Sort by VOG
        sort -k3,3 $TMPDIR/contig_abund.sorted > $TMPDIR/tmp1

        #Sort the virlin by vog
        sort -k1,1 $TMPDIR/viruses > $TMPDIR/tmp2

        #Join virlin and vog_contig_abundance
        join -1 3 -2 1 $TMPDIR/tmp1 $TMPDIR/tmp2 -t $'\t' | sort -k2,2 -u > $TMPDIR/tmp3


        python3 {params.arrange} $TMPDIR/tmp3 > {output.table}

        """
