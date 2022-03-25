#!/bin/bash -eux
 
acc=$1
ref=${2:-/home/reference/NC_045512/NC_045512.2.fa}
 
bam=${acc}.bam
 
# align reads
hisat2 -x ${ref} --no-spliced-alignment --sra-acc $acc --no-unal --threads 1 | samtools view -Sb - | samtools sort - > ${bam}
 
vcf=${acc}.vcf
 
# compute read pileup
bcftools mpileup -a INFO/AD --max-idepth 1000000 --max-depth 1000000 --fasta-ref ${ref} ${bam} \
| bcftools call -o ${vcf} -Ov --ploidy-file >(echo '* * * * 1') --keep-alts --variants-only --multiallelic-caller  
 
# collect average coverage
bedtools genomecov -d -ibam ${bam} | awk 'BEGIN {sum=0}; {sum+=$3}; END{print sum/NR}' < ${bam}.avg_cov