#!/bin/bash

source /com/extra/python/2.7/load.sh

python SNPdensity.py -i Clover_SP_dbSNP_V1.1.clean.all.vcf -o snp_density.csv -b $1 -d 200
