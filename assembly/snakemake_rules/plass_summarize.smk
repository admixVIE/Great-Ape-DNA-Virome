rule summarize_plass_genes:
    input:
        original=expand("../results/{{experiment}}/04b_plass/{sample}.faa", sample = SAMPLES),
        mmseq=expand("../results/{{experiment}}/06c_selected/{sample}_mmseqs.tsv", sample = SAMPLES),
        vogdb=expand("../results/{{experiment}}/06c_selected/{sample}_vog.tsv", sample = SAMPLES),
    output:
        tab="../results/{experiment}/99_summary/plas_genes_summary.tab"
    threads:1
    resources:
        mem_mb=1000,
        time="01:00:00"
    shell:
        """
        tab="    "

        ar1=({input.original})
        ar2=({input.mmseq})
        ar3=({input.vogdb})

        array_length=${{#ar1[@]}}

        for ((i=0; i<$array_length; i++))
        do
        
        #Sample name in this iteration
        sample=$(basename ${{ar1[$i]}})
        
        plass=$(grep -c "^>" ${{ar1[$i]}} || true)
        
        mmseq=$(wc -l ${{ar2[$i]}} | cut -f 1 -d " ")
        
        vog=$(wc -l ${{ar3[$i]}} | cut -f 1 -d " ")
    
        let viral=mmseq+vog || true
        
        #Write results to a file
        echo "${{sample}}${{tab}}${{plass}}${{tab}}Plass" >> {output.tab}
        echo "${{sample}}${{tab}}${{viral}}${{tab}}viral" >> {output.tab}
        done
        """