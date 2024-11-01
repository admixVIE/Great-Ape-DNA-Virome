rule quast:
    input:
        assemblies=expand("../results/{{experiment}}/04_contigs/{sample}.fasta", sample= SAMPLES)
    output:
        stat="../results/{experiment}/04_contigs_stats/report.tsv"
    threads: 4
    resources: mem_mb=2000,
               time="02:00:00"
    shell:
        """
        module load conda
        conda activate quast

        #do the quast run
        quast.py {input.assemblies} -o $TMPDIR -t {threads}

        #Copy the results
        cp $TMPDIR/report.* $(dirname {output.stat})

        #run the multiqc
        conda activate multiqc
        multiqc $TMPDIR/ --outdir $(dirname {output.stat})
        """