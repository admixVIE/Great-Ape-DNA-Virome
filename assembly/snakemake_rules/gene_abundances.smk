rule gene_abundances:
    input:
        genes="../results/{experiment}/05_genes/{sample}.gff",
        fasta="../results/{experiment}/05_genes/{sample}.proteins.fasta",
        bams1="../results/{experiment}/04_bwa/{sample}.bam",
        bais1="../results/{experiment}/04_bwa/{sample}.bam.bai"
    output:
        abund="../results/{experiment}/05a_gene_abundances/{sample}.abundances.tab",
        dictionary="../results/{experiment}/05a_gene_abundances/{sample}.genedict"
    log:
        "log/gene_abundances_{experiment}_{sample}.log"
    threads: 1
    resources: mem_mb=1000,
               time="4-00:00:00"
    shell:
        """
        module load samtools

        #Exit if input files are empty
        if ! [ -s {input.fasta} ] || ! [ -s {input.bams1} ]; then
            touch {output.abund} {output.dictionary}
            exit 0
        fi
        #Get the names of the contigs as in the voglist (after the contigs name, the ordinal number of predicted gene)
        grep "^>" {input.fasta} | cut -f 1 -d " " | cut -d ">" -f 2 > $TMPDIR/contig.names


        #Merge bam files
        #samtools merge $TMPDIR/tmp.bam {input.bams1} input.bams2
        #samtools index $TMPDIR/tmp.bam
        cp {input.bams1} $TMPDIR/tmp.bam
        cp {input.bais1} $TMPDIR/tmp.bam.bai

        #Extract headers with gene names and coordinates
        grep -v "^#" {input.genes} | cut -f 1,4,5  | awk '{{print $1":"$2"-"$3}}' > $TMPDIR/coords.tab

        #Get the gene abundances
        while read p;do

        #skip the length calculation for now in order to count the number of reads mapping to viral genes
        #length=$(echo $p | cut -d ":" -f 2 | tr "-" " " | awk '{{print $2-$1}}')
        #echo -e $p"\t"$(samtools view $TMPDIR/tmp.bam $p | wc -l | awk -v len=$length '{{print $1/len}}') >> $TMPDIR/tmp.abund
        echo -e $p"\t"$(samtools view $TMPDIR/tmp.bam $p | wc -l | awk '{{print $1}}') >> $TMPDIR/tmp.abund
        done < $TMPDIR/coords.tab

        #Adjust the names of the contigs so they match those in voglist
        cut -f 2 $TMPDIR/tmp.abund >> $TMPDIR/abundance.list
        paste -d '\t' $TMPDIR/contig.names $TMPDIR/abundance.list > {output.abund}

        #Create a dictionary of gene names and their coordinates
        paste -d '\t' $TMPDIR/contig.names $TMPDIR/coords.tab > {output.dictionary}
        """
