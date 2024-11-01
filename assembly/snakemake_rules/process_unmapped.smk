rule process_unmapped:
    input:
        bam="../results/{experiment}/04_bwa/{sample}.bam"
    output:
        unmapped="../results/{experiment}/04_unmapped/{sample}.cdhit.fasta.gz",
        clstr="../results/{experiment}/04_unmapped/{sample}.cdhit.clstr.gz"
    threads: 4
    resources:
        mem_mb=4000,
        time="04:00:00"
    shell:
        """
        module load samtools

        #Extract unmapped reads
        samtools view -b -f 12 -F 256 {input.bwa} > $TMPDIR/unmapped.bam
        samtools index $TMPDIR/unmapped.bam

        samtools sort -n $TMPDIR/unmapped.bam -o $TMPDIR/unmapped_nsorted.bam

        samtools fastq -1 $TMPDIR/f1.fastq.gz -2 $TMPDIR/f2.fastq.gz -s $TMPDIR/s.fastq.gz $TMPDIR/unmapped_nsorted.bam

        cat $TMPDIR/*.fastq.gz > $TMPDIR/all.fastq.gz


        #Cluster the reads
        module load cdhit
        cd-hit -i $TMPDIR/all.fastq.gz -o $TMPDIR/{wildcards.sample}.cdhit.fasta -M 4000 -T {threads} -d 0

        #Copy the results
        gzip $TMPDIR/{wildcards.sample}.cdhit.fasta
        gzip $TMPDIR/{wildcards.sample}.cdhit.clstr

        cp $TMPDIR/{wildcards.sample}.cdhit.fasta {output.unmapped}
        cp $TMPDIR/{wildcards.sample}.cdhit.clstr {output.clstr}

        #Clean the tmp dir
        rm -r $TMPDIR

        """
