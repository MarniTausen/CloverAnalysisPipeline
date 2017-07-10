#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
     -T SelectVariants \
     -o /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.MAF.vcf \
     --variant /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.clean.vcf \
     --restrictAllelesTo BIALLELIC -select "AF>0.05"
