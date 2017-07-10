#!/bin/bash

source /com/extra/python/2.7/load.sh

python SFS.py -i Clover_SP_dbSNP_V1.1.clean.vcf -o sfs.csv -d 200 -s 197
