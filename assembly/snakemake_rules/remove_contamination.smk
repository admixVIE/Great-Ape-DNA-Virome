rule remove_contamination:
    input:
        swabs="../results/{experiment}/06_vog_map/swabs.voglist",
        voglist="../results/{experiment}/06_vog_map/{sample}.voglist",
        voglineage="../results/{experiment}/06_vog_map/{sample}.voglin"
    output:
        voglist="../results/{experiment}/06b_vog_map_filtered/{sample}.voglist",
        voglineage="../results/{experiment}/06b_vog_map_filtered/{sample}.voglin",
        viruses="../results/{experiment}/06b_vig_map_filtered/{sample}.virlin"
    threads: 1,
    resources: mem_mb = 1000,
               time = "01:00:00"
    params:
        lca="/lisc/scratch/cube/Viruses_of_great_apes/db/VOGDB2017/vog.lca.tsv"
    shell:
        """
        cut -f 1 {input.swabs} > $TMPDIR/swabs_vogs

        #Remove VOGs that are found in swabs
        grep -f $TMPDIR/swabs_vogs -v {input.voglist} > {output.voglist}

        grep -f $TMPDIR/swabs_vogs -v {input.voglineage} > {output.voglineage}

        

        """
