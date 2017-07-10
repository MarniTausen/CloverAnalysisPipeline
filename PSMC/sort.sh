#!/bin/bash

source /com/extra/samtools/1.3/load.sh
source /com/extra/java/8/load.sh

samtools sort -o $2 $1
