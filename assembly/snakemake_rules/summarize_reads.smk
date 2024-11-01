localrules: summarize_reads
ruleorder: summarize_reads > summarize_reads_SE
rule summarize_reads:
    input:
        raw1=expand("../data/{{experiment}}/{sample}_1.fastq.gz", sample = SAMPLES),
        raw2=expand("../data/{{experiment}}/{sample}_2.fastq.gz", sample = SAMPLES),
        qc1=expand("../results/{{experiment}}/02_trimmed_reads/{sample}_1.fastq.gz", sample = SAMPLES),
        qc2=expand("../results/{{experiment}}/02_trimmed_reads/{sample}_2.fastq.gz", sample = SAMPLES),
        #kra1=expand("../results/{{experiment}}/03_viral_reads/{sample}_1.fastq.gz", sample = SAMPLES),
        #kra2=expand("../results/{{experiment}}/03_viral_reads/{sample}_2.fastq.gz", sample = SAMPLES),
        #vir=expand("../results/{{experiment}}/07_tables/viruses/{sample}.tab", sample = SAMPLES)
    output:
        tab="../results/{experiment}/99_summary/read_summary.tab",
        #plot="../results/{experiment}/99_summary/read_stats.pdf"
    params:
        plotscript="analysis/read_stats.R"
    #log:
    #    "log/summarize_reads_{experiment}.log"
    threads:1
    resources: mem_mb=2000,
               time="1:00:00"
    shell:
        """
        raw1=({input.raw1})
        raw2=({input.raw2})
        qc1=({input.qc1})
        qc2=({input.qc2})

        tab="    "

        array_length=${{#raw1[@]}}

        for ((i=0; i<$array_length; i++))
        do

        #sample name in this iteration
        sample=$(basename ${{raw1[$i]}})
        #echo "1"

        #Number of raw reads
        rawnr1=$(zgrep -c "^+$" ${{raw1[$i]}} || true)
        rawnr2=$(zgrep -c "^+$" ${{raw2[$i]}} || true)
        let "rawnr = rawnr1 + rawnr2"
        #echo "2"

        #Number of QC filtered reads
        qcnr1=$(zgrep -c "^+$" ${{qc1[$i]}} || true)
        qcnr2=$(zgrep -c "^+$" ${{qc2[$i]}} || true)
        let "qcnr = qcnr1 + qcnr2"
        #echo "3"

        #Write to the output file
        echo "${{sample}}${{tab}}${{rawnr}}${{tab}}raw" >> {output.tab}
        echo "${{sample}}${{tab}}${{qcnr}}${{tab}}QC_filter" >> {output.tab}
        #echo "7"
        done

       
        """



rule summarize_reads_SE:
    input:
        raw1=expand("../data/{{experiment}}/{sample}_1.fastq.gz", sample = SAMPLES),
        qc1=expand("../results/{{experiment}}/02_trimmed_reads/{sample}_1.fastq.gz", sample = SAMPLES),
        #kra1=expand("../results/{{experiment}}/03_viral_reads/{sample}_1.fastq.gz", sample = SAMPLES),
        #vir=expand("../results/{{experiment}}/07_tables/viruses/{sample}.tab", sample = SAMPLES)
    output:
        tab="../results/{experiment}/99_summary/read_summary.tab",
        #plot="../results/{experiment}/99_summary/read_stats.pdf"
    params:
        plotscript="analysis/read_stats.R"
    #log:
    #    "log/summarize_reads_{experiment}.log"
    threads:1
    resources: mem_mb=2000,
               time="01:00:00"
    shell:
        """
        module load R

        raw1=({input.raw1})
        qc1=({input.qc1})
        
        tab="    "

        array_length=${{#raw1[@]}}

        for ((i=0; i<$array_length; i++))
        do

        #sample name in this iteration
        sample=$(basename ${{raw1[$i]}})

        #Number of raw reads
        rawdir=$(dirname ${{raw1[$i]}})
        raw2=${{rawdir}}/${{sample}}_2.fastq.gz
        rawnr1=$(zgrep -c "^+$" ${{raw1[$i]}} || true)
        if [ -s $raw2 ]; then
            rawnr2=$(zgrep -c "^+$" ${{raw2}} || true)
        else
            rawnr2=0
        fi
        let "rawnr = rawnr1 + rawnr2"

        #Number of QC filtered reads
        qcdir=$(dirname ${{qc1[$i]}})
        qc2=${{qcdir}}/${{sample}}_2.fastq.gz
        qcnr1=$(zgrep -c "^+$" ${{qc1[$i]}} || true)
        if [ -s $qc2 ]; then
            qcnr2=$(zgrep -c "^+$" ${{qc2}} || true )
        else
            qcnr2=0
        fi
        let "qcnr = qcnr1 + qcnr2"

      
        #Write to the output file
        echo "${{sample}}${{tab}}${{rawnr}}${{tab}}raw" >> {output.tab}
        echo "${{sample}}${{tab}}${{qcnr}}${{tab}}QC_filter" >> {output.tab}

        done


        """
