rule bamtools:
    input:
        bam="../results/{experiment}/04_bwa/{sample}.bam",
    output:
        stat="../results/{experiment}/04_bwa_stats/{sample}.stat",
    threads: 1
    resources: mem_mb=1000,
               time="01:00:00"
    shell:
        """
        module load bamtools
        bamtools stats -in {input.bam} > $TMPDIR/tmp.stat

        cp $TMPDIR/tmp.stat {output.stat}
        rm -r $TMPDIR
        """