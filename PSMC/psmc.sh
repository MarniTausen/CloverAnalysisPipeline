#!/bin/bash

source /com/extra/psmc/2012-11-19/load.sh

fq2psmcfa $1 > ${1%.*}.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${1%.*}.psmc ${1%.*}.psmcfa

