#localrules: trimmomatic

ruleorder: trimmomatic > trimmomatic_SE

rule trimmomatic:
    input:
        R1="../data/{experiment}/{sample}_1.fastq.gz",
        R2="../data/{experiment}/{sample}_2.fastq.gz"
    output:
        G1="../results/{experiment}/02_trimmed_reads/{sample}_1.fastq.gz",
        G2="../results/{experiment}/02_trimmed_reads/{sample}_2.fastq.gz",
        T1="../results/{experiment}/02_trimmed_reads/{sample}_1.trimmed.fastq.gz",
        T2="../results/{experiment}/02_trimmed_reads/{sample}_2.trimmed.fastq.gz"
    log:
        "log/trimmomatic_{experiment}_{sample}.log"
    threads: 2
    resources: mem_mb=2000,
               time="2:00:00"
    shell:
        """
        module load trimmomatic

        cp /lisc/app/trimmomatic/0.39/adapters/TruSeq3-PE-2.fa .


        trimmomatic PE -threads {threads} {input.R1} {input.R2} {output.G1} /dev/null {output.G2} /dev/null ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1:true MINLEN:50
        trimmomatic PE -threads {threads} {input.R1} {input.R2} {output.T1} /dev/null {output.T2} /dev/null ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1:true SLIDINGWINDOW:4:20 LEADING:15 TRAILING:15 MINLEN:50
        """






rule trimmomatic_SE:
    input:
        R1="../data/{experiment}/{sample}_1.fastq.gz",
    output:
        G1="../results/{experiment}/02_trimmed_reads/{sample}_1.fastq.gz",
        T1="../results/{experiment}/02_trimmed_reads/{sample}_1.trimmed.fastq.gz",
    log:
        "log/trimmomatic_SE_{experiment}_{sample}.log"
    threads: 8
    resources: mem_mb=2000,
               time="04:00:00"
    shell:
        """
        module load trimmomatic

        cp /lisc/app/trimmomatic/0.39/adapters/TruSeq3-SE.fa .


        trimmomatic SE -threads {threads} {input.R1} {output.G1}  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1:true MINLEN:50
        trimmomatic SE -threads {threads} {input.R1} {output.T1}  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1:true SLIDINGWINDOW:4:20 LEADING:15 TRAILING:15 MINLEN:50
        """
