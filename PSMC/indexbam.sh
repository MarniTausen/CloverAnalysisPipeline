#!/bin/bash

source /com/extra/samtools/1.3/load.sh

samtools index $1

mv $1.bai $2
