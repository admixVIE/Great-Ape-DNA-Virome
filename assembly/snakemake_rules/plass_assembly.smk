ruleorder: plass > plass_SE

rule plass:
    input:
        R1="../results/{experiment}/02_trimmed_reads/{sample}_1.fastq.gz",
        R2="../results/{experiment}/02_trimmed_reads/{sample}_2.fastq.gz"
    output:
        assembly="../results/{experiment}/04b_plass/{sample}.faa"
    threads: 6
    resources:
        mem_mb=50000,
        time="05:00:00"
    shell:
        """
        module load conda
        conda activate plass

        plass assemble {input.R1} {input.R2} {output.assembly} $TMPDIR --min-length 20 --split-memory-limit 50G
        """

rule plass_SE:
    input:
        R1="../results/{experiment}/02_trimmed_reads/{sample}_1.fastq.gz"
    output:
        assembly="../results/{experiment}/04b_plass/{sample}.faa"
    threads: 6
    resources:
        mem_mb=50000,
        time="05:00:00"
    shell:
        """
        module load conda
        conda activate plass

        plass assemble {input.R1} {output.assembly} $TMPDIR --min-length 20 --split-memory-limit 50G
        """
