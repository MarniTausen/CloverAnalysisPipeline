#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
     -T SelectVariants \
     -o /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.clean.all.vcf \
     --variant /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.final.all.vcf \
     --restrictAllelesTo ALL -select "MQ>30.00 && DP>200 && QUAL>20.00"
