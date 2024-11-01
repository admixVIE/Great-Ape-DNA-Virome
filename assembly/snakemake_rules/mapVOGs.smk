rule hmmscan:
    input:
        fasta="../results/{experiment}/05a_filtered_genes/{sample}.proteins.fasta"
    output:
        hmmscan="../results/{experiment}/06_vog_map/{sample}.tblout",
        voglist="../results/{experiment}/06_vog_map/{sample}.voglist",
        voglineage="../results/{experiment}/06_vog_map/{sample}.voglin",
        viruses="../results/{experiment}/06_vog_map/{sample}.virlin"
    log:
        "log/hmmscan_{experiment}_{sample}.log"
    threads: 1
    resources: mem_mb=16000,
               time="3-00:00:00"
    params:
        #vogdb="/proj/VIROINF/db/VOGdb/vog_database",
        vogdb="/lisc/scratch/cube/Viruses_of_great_apes/db/VOGDB217/vog217",
        extract="bin/extract_vogs_hmmscan.py",
        lca="/lisc/scratch/cube/Viruses_of_great_apes/db/VOGDB217/vog.lca.tsv",
        phageterms="etc/phage_lineage_terms.txt",
        ecutoff=1
    shell:
        """
        module unload python3
        module load hmmer python3/

        #Check if file is not empty
        if [ -s {input.fasta} ]; then

        #Run hmmscan
        hmmsearch --tblout {output.hmmscan} --noali -o /dev/null {params.vogdb} {input.fasta}

        #If there are not matches in the vogdb, touch the output files and exit
        if  [ $(grep -vc "^#" {output.hmmscan}) -eq 0 ]; then
             touch {output.voglist} {output.voglineage} {output.viruses}
             exit 0
        fi

		#Sort by the query sequence and evalue, change order of columns for vog and seq
		grep -v "^#" {output.hmmscan} | sort -b -k1,1 -k5,5g | awk '{{t=$1; $1=$3; $3=t; print}}' > $TMPDIR/tmp.tblout

        #Extract VOGs
        python3 {params.extract} $TMPDIR/tmp.tblout {params.ecutoff} {params.lca} > {output.voglist}

        #Get the lineages of the vogs
        cut -f 1 {output.voglist} > $TMPDIR/vog
        while read -r pattern
        do
        grep "$pattern" {params.lca}
        done < $TMPDIR/vog > {output.voglineage}

        #grep -f $TMPDIR/vog {params.lca} > {output.voglineage}

        #Get the lineages of nonphage viruses
        grep -v -f {params.phageterms} {output.voglineage} | cut -f 4 | tr " " "_" | sort -u > {output.viruses} || true

    else
        touch {output.hmmscan} {output.voglist} {output.voglineage} {output.viruses}
        fi
        #END
        """
