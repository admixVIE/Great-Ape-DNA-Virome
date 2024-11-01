ruleorder: kraken > kraken_SE

rule kraken:
    input:
        R1="../results/{experiment}/02_trimmed_reads/{sample}_1.fastq.gz",
        R2="../results/{experiment}/02_trimmed_reads/{sample}_2.fastq.gz"
    output:
        V1="../results/{experiment}/03_viral_reads/{sample}_1.fastq.gz",
        V2="../results/{experiment}/03_viral_reads/{sample}_2.fastq.gz",
        report="../results/{experiment}/03_kraken_report/{sample}.report.txt"
    params:
        extract_py="bin/extract_kraken_reads.py",
        krakendb="/scratch/cube/Viruses_of_great_apes/db/kraken-all_gorilla"
        #krakendb="/proj/Viruses_of_great_apes/tmp/kraken2_database/kraken2-all_refseq"
    log:
        "log/kraken_{experiment}_{sample}.log"
    threads: 8
    resources: mem_mb=100000,
               time="2-00:00:00"
    shell:
        """
        module unload python3
        module load conda python3/3.10
        conda activate kraken2-2.1.2

        #Run kraken2 analysis
        kraken2 --paired \
                --classified-out $TMPDIR/cseqs#.fq \
                --unclassified-out $TMPDIR/useqs#.fq \
                {input.R1} {input.R2} \
                --threads {threads} \
                --memory-mapping \
                --db {params.krakendb} \
                --fastq-input \
                --gzip-compressed \
                --use-names \
                -report {output.report} \
                > $TMPDIR/kraken.output

        conda deactivate
        module unload conda
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if [ -s $TMPDIR/cseqs_1.fq ]; then

        #Extract reads mapped to viral genomes
        python3 {params.extract_py} -k $TMPDIR/kraken.output \
                                    -s $TMPDIR/cseqs_1.fq \
                                    -s2 $TMPDIR/cseqs_2.fq \
                                    -o $TMPDIR/vseqs_1.fq \
                                    -o2 $TMPDIR/vseqs_2.fq \
                                    --fastq-output \
                                    --include-children \
                                    -r {output.report} \
                                    -t 10239
       fi
       #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #Combine unclassified sequences and viral sequences into a single output file
        touch $TMPDIR/vseqs_1.fq $TMPDIR/useqs_1.fq $TMPDIR/vseqs_2.fq $TMPDIR/useqs_2.fq
        cat $TMPDIR/vseqs_1.fq $TMPDIR/useqs_1.fq | gzip > {output.V1}
        cat $TMPDIR/vseqs_2.fq $TMPDIR/useqs_2.fq | gzip > {output.V2}
        """





rule kraken_SE:
    input:
        R1="../results/{experiment}/02_trimmed_reads/{sample}_1.fastq.gz",
    output:
        V1="../results/{experiment}/03_viral_reads/{sample}_1.fastq.gz",
        report="../results/{experiment}/03_kraken_report/{sample}.report.txt"
    params:
        extract_py="bin/extract_kraken_reads.py",
        krakendb="/scratch/cube/Viruses_of_great_apes/db/kraken-all_gorilla"
        #krakendb="/proj/Viruses_of_great_apes/tmp/kraken2_database/kraken2-all_refseq"
    log:
        "log/kraken_SE_{experiment}_{sample}.log"
    threads: 8
    resources: mem_mb=100000,
               time="2-00:00:00"
    shell:
        """
        module unload python3
        module load conda python3/3.10
        conda activate kraken2-2.1.2


        #Run kraken2 analysis
        kraken2 --classified-out $TMPDIR/cseqs.fq \
                --unclassified-out $TMPDIR/useqs.fq \
                {input.R1}\
                --threads {threads} \
                --memory-mapping \
                --db {params.krakendb} \
                --fastq-input \
                --gzip-compressed \
                --use-names \
                -report {output.report} \
                > $TMPDIR/kraken.output

        conda deactivate
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if [ -s $TMPDIR/cseqs_1.fq ]; then

        #Extract reads mapped to viral genomes
        python3 {params.extract_py} -k $TMPDIR/kraken.output \
                                    -s $TMPDIR/cseqs.fq \
                                    -o $TMPDIR/vseqs.fq \
                                    --fastq-output \
                                    --include-children \
                                    -r {output.report} \
                                    -t 10239
       fi
       #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #Combine unclassified sequences and viral sequences into a single output file
        touch $TMPDIR/vseqs.fq $TMPDIR/useqs.fq
        cat $TMPDIR/vseqs.fq $TMPDIR/useqs.fq | gzip > {output.V1}
        """
