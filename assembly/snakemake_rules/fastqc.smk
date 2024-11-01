localrules: fastQC

rule fastQC:
    input:
        "../data/{experiment}/{sample}_{read}.fastq.gz"
    output:
        "../results/{experiment}/01_fastQC/{sample}_{read}_fastqc.html"
    log:
        "log/fastqc_{experiment}_{sample}_{read}.log"
    threads: 1
    resources: mem_mb=1000,
               time="1:00:00"
    shell:
        """
        module load fastqc

        fastqc {input} -o $(dirname {output})
        """
