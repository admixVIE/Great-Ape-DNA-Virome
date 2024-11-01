rule assemble_SE:
    input:
        R1="../results/{experiment}/02_trimmed_reads/{sample}_1.fastq.gz",
    output:
        contigs="../results/{experiment}/04_contigs/{sample}.fasta",
        bam="../results/{experiment}/04_bwa/{sample}.bam",
        bai="../results/{experiment}/04_bwa/{sample}.bam.bai"
    log:
        "log/assemble_SE_{experiment}_{sample}.log"
    threads: 8
    resources: mem_mb=150000,
               time="5-00:00:00"
    shell:
        """
        module load spades bwa samtools
        which metaviralspades.py
        which metaspades.py


        #Check if the input files are empty. If yes, skip the whole execution
        testvar=$(wc -l {input.R1} | cut -f 1 -d " ")
        if ! [ $testvar == 0 ]; then

        ########################### Assemble by metaviralspades and map reads back ####################

        #Create temporary output directories
        mkdir -p $TMPDIR/mvs $TMPDIR/metaspades $TMPDIR/bwa/reference $TMPDIR/unmapped $TMPDIR/bwa2/reference

        #Assemble with metaviralspades
        spades.py -s {input.R1} -o $TMPDIR/mvs -t {threads} --phred-offset 33

        #Map reads back to the assemblies

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if [ -s $TMPDIR/mvs/contigs.fasta ]; then
        cp $TMPDIR/mvs/contigs.fasta $TMPDIR/bwa/reference/

        bwa index $TMPDIR/bwa/reference/contigs.fasta

        bwa mem $TMPDIR/bwa/reference/contigs.fasta {input.R1} -t {threads} -v 2 \
            | samtools view -b > $TMPDIR/bwa/tmp.bam

        samtools sort $TMPDIR/bwa/tmp.bam -o $TMPDIR/bwa/contigs.bam

        rm $TMPDIR/bwa/tmp.bam

        samtools index $TMPDIR/bwa/contigs.bam
        fi
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



        ######################### Extract unmapped reads ##############################################

        #Extract unmapped reads

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if [ -s $TMPDIR/bwa/contigs.bam ]; then

        samtools view -b -f 4 -F 256 $TMPDIR/bwa/contigs.bam > $TMPDIR/unmapped/unmapped.bam

        samtools index $TMPDIR/unmapped/unmapped.bam

        #Create fastq files
        samtools sort -n $TMPDIR/unmapped/unmapped.bam -o $TMPDIR/unmapped/nsorted.bam
        samtools fastq $TMPDIR/unmapped/nsorted.bam | gzip > $TMPDIR/unmapped/fq1.fastq.gz
        fi
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        ####################   Assemble unmapped read with metaspades and map reads back #############

        #Assemble
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if [ -s $TMPDIR/unmapped/fq1.fastq.gz ]; then
        spades.py -s $TMPDIR/unmapped/fq1.fastq.gz -t {threads} -o $TMPDIR/metaspades --phred-offset 33
      else
        spades.py -s {input.R1} -t {threads} -o $TMPDIR/metaspades --phred-offset 33
        cp {input.R1} $TMPDIR/unmapped/fq1.fastq.gz
        fi
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #Map reads back

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if [ -s $TMPDIR/metaspades/contigs.fasta ]; then
        cp $TMPDIR/metaspades/contigs.fasta $TMPDIR/bwa2/reference/

        bwa index $TMPDIR/bwa2/reference/contigs.fasta

        bwa mem $TMPDIR/bwa2/reference/contigs.fasta $TMPDIR/unmapped/fq1.fastq.gz -t {threads} -v 2 | samtools view -b > $TMPDIR/bwa2/tmp.bam

        samtools sort $TMPDIR/bwa2/tmp.bam -o $TMPDIR/bwa2/contigs.bam

        samtools index $TMPDIR/bwa2/contigs.bam

        rm $TMPDIR/bwa2/tmp.bam
        fi
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



       ######################### Merge output files #####################################

       #Merge contigs

       #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       if [ -s $TMPDIR/mvs/contigs.fasta ] && [ -s $TMPDIR/metaspades/contigs.fasta ]; then
       cat $TMPDIR/mvs/contigs.fasta $TMPDIR/metaspades/contigs.fasta > {output.contigs}
      elif ! [ -s $TMPDIR/mvs/contigs.fasta ] && [ -s $TMPDIR/metaspades/contigs.fasta ]; then
       cat $TMPDIR/metaspades/contigs.fasta > {output.contigs}
      elif [ -s $TMPDIR/mvs/contigs.fasta ] && ! [ -s $TMPDIR/metaspades/contigs.fasta ]; then
       cat $TMPDIR/mvs/contigs.fasta > {output.contigs}
      else
       touch {output.contigs}
       fi
       #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

       #Merge bam files

       #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       if [ -s $TMPDIR/bwa/contigs.bam ] && [ -s $TMPDIR/bwa2/contigs.bam ]; then
       samtools merge {output.bam} $TMPDIR/bwa/contigs.bam $TMPDIR/bwa2/contigs.bam
       samtools index {output.bam}
     elif [ -s $TMPDIR/bwa/contigs.bam ] && ! [ -s $TMPDIR/bwa2/contigs.bam ]; then
       cp $TMPDIR/bwa/contigs.bam  {output.bam}
       samtools index {output.bam}
     elif ! [ -s $TMPDIR/bwa/contigs.bam ] && [ -s $TMPDIR/bwa2/contigs.bam ]; then
        cp $TMPDIR/bwa2/contigs.bam  {output.bam}
        samtools index {output.bam}
     else
       touch {output.bam} {output.bai}
       fi
       #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


       else
        touch {output.bam}
        touch {output.contigs}
       fi

        """
