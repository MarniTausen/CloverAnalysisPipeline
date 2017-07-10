#!/bin/bash

source /com/extra/bwa/0.7.5a/load.sh
source /com/extra/samtools/1.3/load.sh

bwa mem -t 8 -R '@RG\tID:1\tSM:1' /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
    $1 $2 | samtools view -Sb - > $5

bwa mem -t 8 -R '@RG\tID:1\tSM:1' /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
    $3 $4 | samtools view -Sb - > $6
