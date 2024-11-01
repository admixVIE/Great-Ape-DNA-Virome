rule summary:
    input:
        raw1=expand("../data/{{experiment}}/{sample}_1.fastq.gz",sample=SAMPLES),
        raw2=expand("../data/{{experiment}}/{sample}_2.fastq.gz",sample=SAMPLES),
        trim1=expand("../results/{{experiment}}/02_trimmed_reads/{sample}_1.fastq.gz",sample=SAMPLES),
        trim2=expand("../results/{{experiment}}/02_trimmed_reads/{sample}_2.fastq.gz", sample=SAMPLES),
        viral1=expand("../results/{{experiment}}/03_viral_reads/{sample}_1.fastq.gz", sample=SAMPLES),
        viral2=expand("../results/{{experiment}}/03_viral_reads/{sample}_2.fastq.gz", sample=SAMPLES),
        contigs=expand("../results/{{experiment}}/04_contigs/{sample}.fasta",sample=SAMPLES)
    output:
        summary="../results/{experiment}/00_summary/summary.tab"
    log:
        "log/summary_{experiment}.log"
    resources: mem_mb=2000,
               time="1:00:00"
    threads: 1
    shell:
        """
        #Get arrays
        raw1=({input.raw1})
        trim1=({input.trim1})
        viral1=({input.viral1})
        contigs=({input.contigs})

        #Table header
        echo -e "sampleName\trawReads\ttrimmedReads\tviralReads\tcontigs" > {output.summary}

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        for ((i=0; i<${{#raw1[@]}}; i++))
        do

        name=$(basename ${{contigs[$i]}} .fasta)

        nr_raw=$(zgrep -c "^@" ${{raw1[$i]}} || true)

        nr_trim=$(zgrep -c "^@" ${{trim1[$i]}} || true)

        nr_viral=$(zgrep -c "^@" ${{viral1[$i]}} || true)

        nr_contigs=$(grep -c "^>" ${{contigs[$i]}} || true)

        echo -e "$name\t${{nr_raw}}\t${{nr_trim}}\t${{nr_viral}}\t${{nr_contigs}}" >> {output.summary}
        done
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        """
