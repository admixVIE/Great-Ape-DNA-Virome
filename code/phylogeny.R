## SNP calling

module load bcftools
module load conda
conda activate /lisc/user/haemmerle/miniconda3_envs/VIRAL_Python3_CONDA/
  
  for i in $(find /lisc/scratch/admixlab/michelle/sequencing_runs/mappings_DNA_virome_paper/L3859 -name "*_alnSE_sorted_rmdup_corrected.bam"); do
sample=$(echo "${i}" | awk -F".bam" '{print $1}'| awk -F"/" '{print $NF}')
id=$(echo "${i}" | awk -F"/" '{print $NF}' | awk -F"alnSE" '{print $1}')
freebayes -f /lisc/scratch/admixlab/michelle/sequencing_runs/refgenomes/HBV_genomes_Locarnini_paper_2021/AF305327.fasta  --report-monomorphic --min-alternate-count 5 --min-coverage 5 -m 30 -F 0.9 --ploidy 1 ${i} > ${sample}_freebayes.vcf
bcftools filter -O v -o ${sample}_freebayes_qual.vcf -s LOWQUAL -i '(INFO/DP>=10)' ${sample}_freebayes.vcf
bcftools norm -f /lisc/scratch/admixlab/michelle/sequencing_runs/refgenomes/HBV_genomes_Locarnini_paper_2021/AF305327.fasta  -m -snps ${sample}_freebayes_qual.vcf -O v -o ${sample}_freebayes_qual_norm.vcf
awk 'length($4) == length($5) || /^#/ {print $0}' ${sample}_freebayes_qual_norm.vcf > ${sample}_freebayes_qual_norm_noInDel.vcf
bcftools view ${sample}_freebayes_qual_norm_noInDel.vcf -Oz -o ${sample}_freebayes_qual_norm_noInDel.vcf.gz
bcftools index ${sample}_freebayes_qual_norm_noInDel.vcf.gz
bcftools consensus --exclude 'FILTER="LOWQUAL"' -f /lisc/scratch/admixlab/michelle/sequencing_runs/refgenomes/HBV_genomes_Locarnini_paper_2021/AF305327.fasta  ${sample}_freebayes_qual_norm_noInDel.vcf.gz -o ${sample}_freebayes_qual_norm_noInDel_bcftools.fasta
sed -i "s/>.*/>$id/g" ${sample}_freebayes_qual_norm_noInDel_bcftools.fasta
done

## Alignment with MAFFT
module load mafft

mafft --thread 8 --auto   HBV_new_combined.fasta  >  MAFFT_aligned_HBV_241014.fasta 


## IQtree
##Submission Varialbles
date
pwd -L
#module purge
set -o errexit
module load iqtree

iqtree2 -s /lisc/scratch/admixlab/michelle/sequencing_runs/phylogeny_DNA_virome_paper/MAFFT_aligned_HBV_241014.fas  -m MFP -nt 8 --prefix /lisc/scratch/admixlab/michelle/sequencing_runs/phylogeny_DNA_virome_paper/MAFFT_aligned_90PD_HBV_241014.fas -b 1000


