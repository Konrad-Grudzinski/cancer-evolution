#!/bin/bash

INFILE=$1
VCF=$2.vcf

# align
mafft --thread 4 --auto --quiet --add ${INFILE} NC_045512.2.fa | tee temp_aligned.fa > /dev/null

# create SAM
'/mnt/c/Users/Konrad Grudzinski/OneDrive - University of Glasgow/Computing/4th Year/Individual Project/Source/hisat2-2.2.1/hisat2' -x NC_045512.2.fa -f temp_aligned.fa -S temp_all.sam
#hisat2 -x NC_045512.2.fa -f temp_aligned.fa -S temp_all.sam

# sort it
samtools sort temp_all.sam -o temp_all-sorted.bam

# call it
#samtools mpileup -uf NC_045512.2.fa temp_all-sorted.bam | bcftools view - > ${VCF}
#samtools mpileup -uf NC_045512.2.fa temp_all-sorted.bam | '/mnt/c/Users/Konrad Grudzinski/OneDrive - University of Glasgow/Computing/4th Year/Individual Project/Source/bcftools-1.14/bcftools' view - > ${VCF}
samtools mpileup -uf NC_045512.2.fa temp_all-sorted.bam | tee temp_pileup > /dev/null
'/mnt/c/Users/Konrad Grudzinski/OneDrive - University of Glasgow/Computing/4th Year/Individual Project/Source/bcftools-1.14/bcftools' view temp_pileup | tee ${VCF} > /dev/null