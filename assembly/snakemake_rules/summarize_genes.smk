localrules: summarize_genes

rule summarize_genes:
    input:
        original_genes=expand("../results/{{experiment}}/05_genes/{sample}.proteins.fasta", sample=SAMPLES),
        diamond_genes=expand("../results/{{experiment}}/05a_filtered_genes/{sample}.proteins.fasta", sample=SAMPLES),
        vogdb_genes=expand("../results/{{experiment}}/06_vog_map/{sample}.voglist",sample=SAMPLES)
    output:
        tab="../results/{experiment}/99_summary/gene_summary.tab",
        plot="../results/{experiment}/99_summary/gene_stats.pdf"
    params:
        plotscript="analysis/gene_stats.R"
    log:
        "log/summarize_genes_{experiment}.log"
    threads:1
    resources: mem_mb=1000,
               time="1:00:00"
    shell:
        """
        tab="    "
        ar1=({input.original_genes})
        ar2=({input.diamond_genes})
        ar3=({input.vogdb_genes})

        array_length=${{#ar1[@]}}

        for ((i=0; i<$array_length; i++))
        do

        #sample name in this iteration
        sample=$(basename ${{ar1[$i]}})

        #Number of genes predicted by prodigal
        prod=$(grep -c "^>" ${{ar1[$i]}} || true)

        #Number of genes after DIAMOND filtering
        diam=$(grep -c "^>" ${{ar2[$i]}} || true)

        #Number of genes mapped to VOGs
        vog=$(wc -l ${{ar3[$i]}} | cut -f 1 -d " ")

        #Write results to a file
        echo  "${{sample}}${{tab}}${{prod}}${{tab}}prodigal" >> {output.tab}
        echo  "${{sample}}${{tab}}${{diam}}${{tab}}DIAMOND" >> {output.tab}
        echo  "${{sample}}${{tab}}${{vog}}${{tab}}VOG" >> {output.tab}

        done

        #Create plots
        module load R
        Rscript --vanilla {params.plotscript}

        """
