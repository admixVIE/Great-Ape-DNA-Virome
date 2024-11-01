rule prodigal:
    input:
        #contig="../results/{experiment}/04_contigs/{sample}.fasta",
        contig="../results/{experiment}/046_filtered_contigs/{sample}.fasta",
    output:
        fasta="../results/{experiment}/05_genes/{sample}.proteins.fasta",
        nfasta="../results/{experiment}/05_genes/{sample}.genes.fasta",
        genes="../results/{experiment}/05_genes/{sample}.genes"
    log:
        "log/prodigal_{experiment}_{sample}.log"
    threads: 1
    resources: mem_mb=1000,
               time="5:00:00"
    shell:
        """
        module load prodigal

        #Check if file is empty
        if [ -s {input.contig} ]; then
        prodigal -i {input.contig} -o {output.genes} -a {output.fasta} -d {output.nfasta} -p meta

        else
        touch {output.genes} {output.fasta} {output.nfasta}
        fi
        """
