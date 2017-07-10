#!/bin/bash

source /com/extra/vcftools/0.1.14/load.sh

/com/extra/vcftools/0.1.14/bin/vcftools --vcf Clover_SP_dbSNP_V1.1.MAF.vcf --geno-r2 --ld-window-bp 10000 --out LD_10k_window
