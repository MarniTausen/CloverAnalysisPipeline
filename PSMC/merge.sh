#!/bin/bash

source /com/extra/samtools/1.3/load.sh

samtools merge $3 $1 $2
