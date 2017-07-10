#!/bin/bash

source /com/extra/Anaconda-Python/2.2.0-2.7/load.sh
source activate cloverevolution

python2.7 VCFtoGRM.py -i Clover_SP_dbSNP_V1.1.MAF.vcf -o GRM.csv -n 20
