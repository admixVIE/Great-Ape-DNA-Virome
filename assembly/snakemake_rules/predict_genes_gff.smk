rule prodigal_gff:
    input:
        #contig="../results/{experiment}/04_contigs/{sample}.fasta",
        contig="../results/{experiment}/046_filtered_contigs/{sample}.fasta",
    output:
        gff="../results/{experiment}/05_genes/{sample}.gff"
    log:
        "log/prodigal_gff_{experiment}_{sample}.log"
    threads: 1
    resources: mem_mb=1000,
               time="5:00:00"
    shell:
        """
        module load prodigal

        #Check if file is empty
        if [ -s {input.contig} ]; then
        prodigal -i {input.contig} -o {output.gff} -f gff -p meta

        else
        touch {output.gff}
        fi
        """
