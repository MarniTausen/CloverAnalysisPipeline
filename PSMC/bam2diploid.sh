#!/bin/bash

source /com/extra/samtools/1.4.1/load.sh
source /com/extra/bcftools/1.4.1/load.sh

samtools mpileup -C 50 -uf /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta $1 | bcftools call -c | vcfutils.pl vcf2fq -d 5 -D 100 - | gzip > ${1%.*}.fq.gz


