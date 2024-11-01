ruleorder: assemble > assemble_SE

rule assemble:
    input:
        R1="../results/{experiment}/02_trimmed_reads/{sample}_1.fastq.gz",
        R2="../results/{experiment}/02_trimmed_reads/{sample}_2.fastq.gz",
    output:
        contigs="../results/{experiment}/04_contigs/{sample}.fasta",
        bam="../results/{experiment}/04_bwa/{sample}.bam",
        bai="../results/{experiment}/04_bwa/{sample}.bam.bai"
    log:
        "log/assemble_{experiment}_{sample}.log"
    threads: 8
    resources: mem_mb=200000,
               time="5-00:00:00"
    shell:
        """
        module load spades bwa samtools seqtk
        which metaviralspades.py
        which metaspades.py


        #Check if the input files are empty. If yes, skip the whole execution
        testvar=$(wc -l {input.R1} | cut -f 1 -d " ")
        if ! [ $testvar == 0 ]; then

        ########################### Assemble by metaviralspades and map reads back ####################

        #Create temporary output directories
        mkdir -p $TMPDIR/mvs $TMPDIR/metaspades $TMPDIR/bwa/reference $TMPDIR/unmapped $TMPDIR/bwa2/reference

        #Assemble with metaviralspades
        metaviralspades.py -1 {input.R1} -2 {input.R2} -o $TMPDIR/mvs -t {threads} --phred-offset 33 -m 200

        #Get the error corrected reads
        r1=($(ls $TMPDIR/mvs/corrected/*_1.fastq*.fastq.gz))
        r2=($(ls $TMPDIR/mvs/corrected/*_2.fastq*.*.fastq.gz))
        ru=($(ls $TMPDIR/mvs/corrected/*__unpaired*.fastq.gz))
        #Merge error corrected reads
        seqtk mergepe $r1 $r2 | gzip > $TMPDIR/mvs/corrected_reads.fastq.gz
        cat $ru >> $TMPDIR/mvs/corrected_reads.fastq.gz

        #Map reads back to the assemblies

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if [ -s $TMPDIR/mvs/contigs.fasta ]; then
        cp $TMPDIR/mvs/contigs.fasta $TMPDIR/bwa/reference/

        bwa index $TMPDIR/bwa/reference/contigs.fasta

        bwa mem $TMPDIR/bwa/reference/contigs.fasta -p $TMPDIR/mvs/corrected_reads.fastq.gz -t {threads} -v 2 \
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

        samtools view -b -f 12 -F 256 $TMPDIR/bwa/contigs.bam > $TMPDIR/unmapped/unmapped.bam

        samtools index $TMPDIR/unmapped/unmapped.bam

        #Create fastq files
        samtools sort -n $TMPDIR/unmapped/unmapped.bam -o $TMPDIR/unmapped/nsorted.bam
        samtools fastq -1 $TMPDIR/unmapped/fq1.fastq.gz -2 $TMPDIR/unmapped/fq2.fastq.gz -s $TMPDIR/unmapped/singletons.fastq.gz $TMPDIR/unmapped/nsorted.bam
        fi
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        ####################   Assemble unmapped read with metaspades and map reads back #############

        #Assemble
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if [ -s $TMPDIR/unmapped/fq1.fastq.gz ]; then
            if ! [ $(zcat $TMPDIR/unmapped/singletons.fastq.gz | wc -l) == 0 ]; then
                metaspades.py -1 $TMPDIR/unmapped/fq1.fastq.gz -2 $TMPDIR/unmapped/fq2.fastq.gz -s $TMPDIR/unmapped/singletons.fastq.gz -t {threads} -o $TMPDIR/metaspades --phred-offset 33 -m 200 --only-assembler
            else
                metaspades.py -1 $TMPDIR/unmapped/fq1.fastq.gz -2 $TMPDIR/unmapped/fq2.fastq.gz -t {threads} -o $TMPDIR/metaspades --phred-offset 33 -m 200 --only-assembler
            fi
      else
        metaspades.py -1 {input.R1} -2 {input.R2} -t {threads} -o $TMPDIR/metaspades --phred-offset 33
        cp {input.R1} $TMPDIR/unmapped/fq1.fastq.gz
        cp {input.R2} $TMPDIR/unmapped/fq2.fastq.gz
        fi
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #Map reads back

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if [ -s $TMPDIR/metaspades/contigs.fasta ]; then
        cp $TMPDIR/metaspades/contigs.fasta $TMPDIR/bwa2/reference/

        bwa index $TMPDIR/bwa2/reference/contigs.fasta

        bwa mem $TMPDIR/bwa2/reference/contigs.fasta $TMPDIR/unmapped/fq1.fastq.gz $TMPDIR/unmapped/fq2.fastq.gz -t {threads} -v 2 | samtools view -b > $TMPDIR/bwa2/tmp.bam

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
