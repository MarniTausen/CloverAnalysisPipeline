Clover genotype analysis supplementary
================
26/03/2018 - 13:40:57

-   [Introduction](#introduction)
-   [GBS genotyping](#gbs-genotyping)
    -   [Unifiedgenotyper](#unifiedgenotyper)
    -   [Selectvariants](#selectvariants)
    -   [SNP density](#snp-density)
    -   [Site Frequency Spectrum](#site-frequency-spectrum)
    -   [Parse VCF](#parse-vcf)
        -   [LD analysis](#ld-analysis)
    -   [Genetic Relationship Matrix](#genetic-relationship-matrix)
-   [GBS data analysis](#gbs-data-analysis)
    -   [SNP density](#snp-density-1)
    -   [Site Frequency Spectrum (SFS)](#site-frequency-spectrum-sfs)
    -   [Genetic Relationship Matrix (GRM)](#genetic-relationship-matrix-grm)
    -   [Heterozygosity](#heterozygosity)
-   [PSMC Pipeline](#psmc-pipeline)
    -   [Simulations](#simulations)
-   [Divergent genes pipeline](#divergent-genes-pipeline)
    -   [FASTafilter](#fastafilter)
    -   [dNdS script](#dnds-script)
    -   [Tabixsearch](#tabixsearch)
    -   [vcf2fasta.py](#vcf2fasta.py)

Introduction
============

This document has all of the scripts used for the Clover genotype analysis. The results are included in a separate document. Each step of the process has its own section and code in the correct order it was run.

GBS genotyping
==============

The GBS genotyping of the 200 clover individuals. The main workflow script is a gwf pipeline (<http://gwf.readthedocs.io/en/latest/>), which keeps track of all job dependencies and submits the scripts in the correct order.

``` python
from gwf import Workflow

gwf = Workflow()

## Genotyping

## SNPS
gwf.target("Genotyper", outputs = ['Clover_SP_dbSNP_V1.1.final.vcf'],
           cores=2, memory='32g', walltime="240:00:00", account="NChain") << """
./Unifiedgenotyper.sh
"""

## ALL_SITES
gwf.target("Genotyper_all", outputs = ['Clover_SP_dbSNP_V1.1.final.all.vcf'],
           cores=2, memory='32g', walltime="240:00:00", account="NChain") << """
./Unifiedgenotyperall.sh
"""

## SNPS
gwf.target("SelectVariants",
           inputs=['Clover_SP_dbSNP_V1.1.final.vcf'],
           outputs = ['Clover_SP_dbSNP_V1.1.clean.vcf'],
           cores=1, memory='32g', walltime="08:00:00", account="NChain") << """
./Selectvariants.sh
"""

gwf.target("SelectVariantsMAF5",
           inputs=['Clover_SP_dbSNP_V1.1.clean.vcf'],
           outputs=['Clover_SP_dbSNP_V1.1.MAF.vcf'],
           cores=1, memory='32g', walltime="08:00:00", account="NChain") << """
./SelectvariantsMAF5.sh
"""

## ALL_SITES
gwf.target("SelectVariants_all",
           inputs=['Clover_SP_dbSNP_V1.1.final.all.vcf'],
           outputs = ['Clover_SP_dbSNP_V1.1.clean.all.vcf'],
           cores=1, memory='32g', walltime="08:00:00", account="NChain") << """
./Selectvariantsall.sh
"""

gwf.target("SelectVariantsMAF5_all",
           inputs=['Clover_SP_dbSNP_V1.1.clean.all.vcf'],
           outputs=['Clover_SP_dbSNP_V1.1.MAF.all.vcf'],
           cores=1, memory='32g', walltime="08:00:00", account="NChain") << """
./SelectvariantsMAF5all.sh
"""

## SNP density

gwf.target("SNPdensity", inputs=['Clover_SP_dbSNP_V1.1.clean.all.vcf'],
           outputs = ['snp_density.csv'],
           cores=1, memory='90g', walltime="12:00:00", account="NChain") << """
./SNPdensity.sh 100
"""

## Site Frequency Spectrum

gwf.target("SFS", inputs=['Clover_SP_dbSNP_V1.1.clean.vcf'],
           outputs = ['sfs.csv'],
           cores=1, memory='8g', walltime="02:00:00", account="NChain") << """
./SFS.sh
"""

## Linkage disequilibrium

gwf.target("LD_vcftools", inputs=['Clover_SP_dbSNP_V1.1.MAF.vcf'],
           outputs = ['LD_10k_window.geno.ld'],
           cores=1, memory='4g', walltime="04:00:00", account="NChain") << """
./LD.sh
"""

gwf.target("LD_figures", inputs=['LD_10k_window.geno.ld'],
           outputs=['chr1.png','chr2.png','chr3.png','chr4.png','chr5.png',
                    'chr6.png','chr7.png','chr8.png','chr9.png','chr10.png',
                    'chr11.png','chr12.png','chr13.png','chr14.png','chr15.png',
                    'chr16.png'],
           cores=1, memory='4g', walltime="04:00:00", account="NChain") << """
./LDR.sh
"""

## Genetic Relationship Matrix

gwf.target("GRM", inputs=['Clover_SP_dbSNP_V1.1.MAF.vcf'],
           outputs=['GRM.csv'], memory='2g', account="NChain") << """
./GRM.sh
"""
```

### Unifiedgenotyper

./Unifiedgenotyper.sh, Produces the unfiltered vcf file for only variant sites.

``` python
#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -T UnifiedGenotyper \
     -I /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/bam.list \
     -o /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.final.vcf \
     -stand_call_conf 50 \
     -stand_emit_conf 20.0 \
     --sample_ploidy 2 \
     -nct 12 --genotyping_mode DISCOVERY \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
     --output_mode EMIT_VARIANTS_ONLY \
     -rf BadCigar 2>&1 > log
```

#### Unifiedgenotyper all

./Unifiedgenotyperall.sh, Produces the unfiltered gvcf file which includes all confident sites.

``` python
#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -T UnifiedGenotyper \
     -I /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/bam.list \
     -o /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.final.all.vcf \
     -stand_call_conf 50 \
     -stand_emit_conf 20.0 \
     --sample_ploidy 2 \
     -nct 12 --genotyping_mode DISCOVERY \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
     --output_mode EMIT_ALL_CONFIDENT_SITES \
     -rf BadCigar 2>&1 > log
```

### Selectvariants

./Selectvariants.sh, Filters all sites which have a Mapping quality less than 30, and a Depth less than 200 and Quality less than 20.

``` python
#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
     -T SelectVariants \
     -o /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.clean.vcf \
     --variant /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.final.vcf \
     --restrictAllelesTo ALL -select "MQ>30.00 && DP>200 && QUAL>20.00"
```

#### Select variants Minor Allele Frequency

./SelectvariantsMAF5.sh, Uses the filtered vcf file and filters for minor allele frequency above 5%

``` python
#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
     -T SelectVariants \
     -o /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.MAF.vcf \
     --variant /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.clean.vcf \
     --restrictAllelesTo BIALLELIC -select "AF>0.05"
```

#### Select variants All sites

./Selectvariantsall.sh, Filters all sites which have a Mapping quality less than 30, and a Depth less than 200 and Quality less than 20 in the gvcf file with all sites.

``` python
#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
     -T SelectVariants \
     -o /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.clean.all.vcf \
     --variant /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.final.all.vcf \
     --restrictAllelesTo ALL -select "MQ>30.00 && DP>200 && QUAL>20.00"
```

#### Select variants Minor Allele Frequency All sites

./SelectvariantsMAF5all.sh, Uses the filtered vcf file and filters for minor allele frequency above 5% in the filtered gvcf file.

``` python
#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
     -T SelectVariants \
     -o /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.MAF.all.vcf \
     --variant /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/WORKFLOW/Clover_SP_dbSNP_V1.1.clean.all.vcf \
     --restrictAllelesTo BIALLELIC -select "AF>0.05"
```

### SNP density

SNP density calculations using the filtered genome vcf sites. Dividing the genome into bins and counting variants and total number of sites. The SNP density is then variants/total number of sites.

``` python
#!/bin/bash

source /com/extra/python/2.7/load.sh

python SNPdensity.py -i Clover_SP_dbSNP_V1.1.clean.all.vcf -o snp_density.csv -b $1 -d 200
```

SNPdensity.py script

``` python
from parseVCF import readVCF
from sys import argv, exit
from optparse import OptionParser

def collectReferenceLengths(header):
    contigs = header.split("contig")[1:]
    contigs = [contig.split(">")[0] for contig in contigs]

    lengths = {}
    
    for contig in contigs:
        ID, length = contig.split(",")
        ID = ID.split("=")[-1]
        length = length.split("=")[-1]
        
        lengths[ID] = int(length)
    
    return lengths

def makeBins(referenceSizes, binsize):

    bins = {}
    
    for chromosome, length in referenceSizes.items():
        bins[chromosome] = [{'sites':0, 'SNP':0, 'DP':0, 'AF':0}
                            for _ in range(length/binsize+1)]

    return bins

def countSNPS(data, bins, binsize, depth):

    for i in range(len(data['CHROM'])):
        chrom, pos, alt, dp, af = data['CHROM'][i], data['POS'][i], data['ALT'][i], data['INFO'][i]['DP'], data['INFO'][i].get('AF', 0)

        bpos = pos/binsize

        if dp<depth: continue

        bins[chrom][bpos]['sites'] += 1
        bins[chrom][bpos]['DP'] += dp
        if alt!=".": 
            bins[chrom][bpos]['SNP'] += 1
            try:
                bins[chrom][bpos]['AF'] += af
            except TypeError:
                print "AF was a string:", af

def SNPdensity(bins):
    SNPdens = {}

    for chromosome, counts in bins.items():
        SNPdens[chromosome] = []
        for count in counts:
            if count['sites']==0:
                SNPdens[chromosome].append('NA')
                continue
            SNPdens[chromosome].append(count['SNP']/float(count['sites']))

    return SNPdens

def printSNPdens(SNPdens, binsize):

    print "%20s %20s %20s" % ('Chromosome', 'Location', 'SNP density')
    
    for chromosome, densities in SNPdens.items():
        i = 0
        for density in densities:
            print "%20s %20s %20s" % (chromosome,
                                       str(i*binsize)+"-"+str((i+1)*binsize),
                                       str(density))
            i += 1

def writeCSV(filename, SNPdens, binsize, bins):

    f = open(filename, 'w')
    
    f.write("%s, %s, %s, %s, %s, %s, %s\n" % ('Chromosome', 'Location', 'SNP density', 'Sites', 'SNPs', 'DP', 'AF'))
    
    for chromosome, densities in SNPdens.items():
        for i, density in enumerate(densities):
            sites = bins[chromosome][i]['sites']
            SNP = bins[chromosome][i]['SNP']
            if sites!=0:
                DP = bins[chromosome][i]['DP']/sites
            else: DP = 0
            if SNP!=0:
                AF = bins[chromosome][i]['AF']/float(SNP)
            else: AF = 0
            f.write("%s, %s, %s, %s, %s, %s, %s\n" % (chromosome,
                                                      str(i*binsize)+"-"+str((i+1)*binsize),
                                                      str(density), str(sites),str(SNP),
                                                      str(DP), str(AF)))

    f.close()

def help():
    print "====== SNP density analysis of a gVCF file ====="
    print "Reads in a gvcf file, strips it down and counts"
    print "the occurrences of SNPS versus the number sites sequenced"
    print "-i <filename.vcf>                 The input vcf file"
    print "-o <filename.csv>                 The output csv file"
    print "-b binsize                        The binsize to be calculated from"
    print "-d read depth                     Filter sites less that the minimum read depth"
    print
    exit()

if __name__=="__main__":

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename.csv>")
    parser.add_option('-i', type="string", nargs=1, dest="input", help="<filename.vcf>")
    parser.add_option('-b', type="int", nargs=1, dest="binsize", help="binsize in kb")
    parser.add_option('-d', type="int", nargs=1, dest="depth", help="read depth filter")
    parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
    options, args = parser.parse_args()

    if options.help!=None:
        help()
    if options.input!=None:
        infile = options.input
    else:
        raise "No input file, please include a vcf file -i <filename.vcf>"

    if options.output!=None:
        outfile = options.output
    else:
        outfile = "output.csv"

    if options.binsize!=None:
        binsize = options.binsize*1000
    else:
        raise "Please provide a binsize: -b binsize"

    if options.depth!=None:
        depth = options.depth
    else:
        print "Warning: Was not provided with a minimum depth filter. Therefore it will not filter any sites"
        print "Please use: -d minimum depth, to filter any sites with a low coverage"
        depth = 0

    print "READING DATA:"
    data, header = readVCF(infile, 9)
    print "RETRIEVING REFERENCE LENGTHS"
    reference = collectReferenceLengths(header)
    print "MAKING BINS"
    bins = makeBins(reference, binsize)
    print "COUNTING SNPS"
    countSNPS(data, bins, binsize, depth)
    print "CALCULATING SNP DENSITIES"
    SNPdens = SNPdensity(bins)
    print "WRITING TO CSV"
    writeCSV(outfile, SNPdens, binsize, bins)
```

Help text from the SNPdensity.py script.

``` python
====== SNP density analysis of a gVCF file =====
Reads in a gvcf file, strips it down and counts
the occurrences of SNPS versus the number sites sequenced
-i <filename.vcf>                 The input vcf file
-o <filename.csv>                 The output csv file
-b binsize                        The binsize to be calculated from
-d read depth                     Filter sites less that the minimum read depth
```

### Site Frequency Spectrum

``` python
#!/bin/bash

source /com/extra/python/2.7/load.sh

python SFS.py -i Clover_SP_dbSNP_V1.1.clean.vcf -o sfs.csv -d 200 -s 197
```

SFS.py script

``` python
from sys import argv
from parseVCF import readVCF
from optparse import OptionParser

def makeSFS(data, depth, samples):

    sfs = [0 for i in range(0, samples*2)]

    for info in data:
        if info['DP']<depth: continue
        if info['AN']<samples: continue
        try:
            sfs[int(info['AC'])] += 1
        except:
            pass

    return sfs

def writeSFS(SFS, outfile):

    f = open(outfile, "w")
    f.write("AN, Count\n")
    for i, c in enumerate(SFS):
        f.write("%i, %i\n" % (i, c))
    f.close()

def help():
    print "====== S Site Frequency spectrum alysis of vcf file ====="
    print "Reads in a vcf file, strips it down and counts"
    print "the occurrences of SNPS versus the number sites sequenced"
    print "-i <filename.vcf>                 The input vcf file"
    print "-o <filename.csv>                 The output csv file"
    print "-d depth                          Read depth to filter by"
    print "-s samples                        Number of samples"

if __name__=="__main__":

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename.csv>")
    parser.add_option('-i', type="string", nargs=1, dest="input", help="<filename.vcf>")
    parser.add_option('-d', type="int", nargs=1, dest="depth", help="read depth filter")
    parser.add_option('-s', type="int", nargs=1, dest="samples", help="number of samples")
    parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
    options, args = parser.parse_args()

    if options.help!=None:
        help()
    if options.input!=None:
        infile = options.input
    else:
        raise "No input file, please include a vcf file -i <filename.vcf>"

    if options.output!=None:
        outfile = options.output
    else:
        outfile = "output.csv"

    if options.depth!=None:
        depth = options.depth
    else:
        print "Warning: Was not provided with a minimum depth filter. Therefore it will not filter any sites"
        print "Please use: -d minimum depth, to filter any sites with a low coverage"
        depth = 0

    if options.samples!=None:
        samples = options.samples
    else:
        raise "Sample size unknown. Please use: -s samples, to filter sites by the number of samples."

    INFO = readVCF(infile, 9)[0]['INFO']
    sfs = makeSFS(INFO, depth, samples)
    writeSFS(sfs, outfile)
```

Parse VCF
---------

parseVCF.py script, which is used by both the SNP density and SFS scripts. Reads and loads a vcf file.

``` python
from sys import argv

def simplify_name(name):
    if "/" in name:
        return name.split("/")[-1].split(".")[0]
    return name

def parseINFO(info):
    stats = {}
    for statistic in info.split(";")[:6]:
        name, value = statistic.split("=")
        try:
            stats[name] = float(value)
        except ValueError:
            stats[name] = value
    return stats

def readVCF(filename, nc=6):
    ''' 
    Reads in a VCF file while discarding most of the data.
    Only reads in, Chromosome (CHROM), position (POS), reference nucleotide (REF),
    alternative nucleotide (ALT), and mapping quality (QUAL). 
    '''

    f = open(filename)
    header = ""
    names = []
    data = {}
    options = {'CHROM': lambda x: data['CHROM'].append(x),
               'POS': lambda x: data['POS'].append(int(x)),
               'ID': lambda x: 1,
               'REF': lambda x: data['REF'].append(x),
               'ALT': lambda x: data['ALT'].append(x),
               'QUAL': lambda x: data['QUAL'].append(float(x)),
               'INFO': lambda x: data['INFO'].append(parseINFO(x)),
               'FORMAT': lambda x: 1,
               'FILTER': lambda x: 1
    }
    for line in f:
        if line[:2]=="##":
            header += line
            continue
        if line[:1]=="#":
            names = line.split("\t")[:nc]
            names = [simplify_name(name) for name in names]
            names[0] = names[0][1:]
            for name in names:
                data[name] = []
            print names
            continue
        for i, point in enumerate(line.split("\t")[:nc]):
            name = names[i]
            if name not in options:
                data[name].append(point)
            else: 
                options[name](point)

    del data['ID']
    
    return data, header


if __name__=="__main__":
    data, header = readVCF(argv[1], int(argv[2]))
    print data
```

### LD analysis

LD analysis, using vcftools with a window size of 10000 bp

``` python
#!/bin/bash

source /com/extra/vcftools/0.1.14/load.sh

/com/extra/vcftools/0.1.14/bin/vcftools --vcf Clover_SP_dbSNP_V1.1.MAF.vcf --geno-r2 --ld-window-bp 10000 --out LD_10k_window
```

./LDR.sh Script loading the Rscript which produces the figures.

``` python
#!/bin/bash

source /com/extra/R/3.3/load.sh

Rscript LD.R
```

LD R script.

``` r
library(ggplot2)
library(dplyr)

LD <- read.table("LD_10k_window.geno.ld", header=TRUE, sep="\t")
#LD <- read.table("TestData.txt", header=TRUE, sep="\t")

#LD$Distance <- abs(LD$POS1-LD$POS2)

LD %>% filter(!is.nan(R.2) & R.2>0 & N_INDV>180) %>%
    mutate(Distance = abs(POS1-POS2)) -> LD

chromosomes <- sort(as.vector(unique(LD$CHR)))

for(chromosome in chromosomes[-1]){
    chrsub <- LD[LD$CHR==chromosome,]
    expmodel <- lm(log(R.2) ~ Distance, data=chrsub)
    r2 <- summary(expmodel)$r.squared
    xseq <- seq(0, 10000, by=0.01)
    reg_fit <- exp(predict(expmodel, list(Distance=xseq)))
    names(reg_fit) <- NULL
    expfit <- data.frame(Distance=xseq, R.2=reg_fit)
    cat("Length of fit:", length(reg_fit), "Length of sequence", length(xseq))
    fig <- ggplot(chrsub, aes(x=Distance, y=R.2)) + geom_point(size=0.3) +
        theme_bw() + labs(title=chromosome, x = "Distance", y = "R-squared linkage") +
        theme(text = element_text(size=8)) +
        annotate("text", label=paste0("R^2 == ", r2), x=5000, y=0.90, parse=TRUE) +
        geom_line(data=expfit, aes(x=Distance, y=R.2), colour="red") +
        scale_x_log10(limits=c(1, 10000))
    ggsave(plot=fig, file=paste0(chromosome, ".png"), width=10, height=4, unit="in")
}
```

Genetic Relationship Matrix
---------------------------

``` python
#!/bin/bash

source /com/extra/Anaconda-Python/2.2.0-2.7/load.sh
source activate cloverevolution

python2.7 VCFtoGRM.py -i Clover_SP_dbSNP_V1.1.MAF.vcf -o GRM.csv -n 20
```

GRM script

``` python
from sys import argv
import numpy as np
import scipy.stats as stats
from optparse import OptionParser

def collectName(string):
    return string.split("/")[-1].split(".")[0]

def getMarker(allele):
    markers = {"0/0": 1.0, "0/1": 0.0, "1/1": -1.0, "./.": None}
    return markers.get(allele.split(":")[0], None)

def readVCFmatrix(filename, n, numpy=True):

    f = open(filename)
    header = ""
    names = []
    matrix = []
    for line in f:
        if line[:2]=="##":
            header += line
            continue
        if line[:1]=="#":
            names = [collectName(name) for name in line.split("\t")[9:]]
            print "Sequence names:", names
            matrix = [[] for i in range(len(names))]
            N = len(names)
            continue
        allele_list = []
        for allele in line.split("\t")[9:]:
            allele_list.append(getMarker(allele))

        if len(filter(lambda x: x!=None, allele_list))>=(N-n):
            for i, allele in enumerate(allele_list):
                if allele!=None:
                    matrix[i].append(allele)
                else:
                    matrix[i].append(0.0)

    if numpy==True:
        matrix = np.array(matrix)

    return matrix, names, header

def replace_column_nans_by_mean(matrix):
    # Set the value of gaps/dashes in each column to be the average of the other values in the column.
    nan_indices = np.where(np.isnan(matrix))
    # Note: bn.nanmean() instead of np.nanmean() because it is a lot(!) faster.
    column_nanmeans = np.nanmean(matrix, axis=0)
    matrix[nan_indices] = np.take(column_nanmeans, nan_indices[1])
    return matrix

def calcGRM(SNP):
    N, M = SNP.shape
    NORM = (SNP-np.mean(SNP, 0))/(np.std(SNP, 0)+0.000000000000001)
    return np.dot(NORM, NORM.T)/M

def help():
    print "====== VCF to GRM (Genetic Relationship Matrix) ====="
    print "Reads in a vcf file as an individual by Marker matrix"
    print "Calculates a GRM using this Marker matrix."
    print "-i <filename.vcf>                 The input vcf file"
    print "-o <filename.csv>                 The output GRM as a csv file"
    print "-n number of missing              How many missing alleles are allowed"

if __name__=="__main__":

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename.csv>")
    parser.add_option('-i', type="string", nargs=1, dest="input", help="<filename.vcf>")
    parser.add_option('-n', type="int", nargs=1, dest="missing", help="number of missing allowed")
    parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
    options, args = parser.parse_args()

    if options.help!=None:
        help()
    if options.input!=None:
        infile = options.input
    else:
        raise "No input file, please include a vcf file -i <filename.vcf>"

    if options.output!=None:
        outfile = options.output
    else:
        outfile = "GRM.csv"

    if options.missing!=None:
        nmissing = int(options.missing)
    else:
        nmissing = 0

    print "READING VCF FILE"
    matrix, names, header = readVCFmatrix(infile, nmissing)

    print "STARTING CALCULATIONS"
    print "Number of individuals", matrix.shape[0]
    M = matrix.shape[1]
    print "Number of Markers", M
    matrix = replace_column_nans_by_mean(matrix)
    GRM = calcGRM(matrix)

    print "SAVING FILE"
    np.savetxt(outfile, GRM, delimiter=", ")
```

GBS data analysis
=================

This section contains all of the R scripts which visualize the results from the GBS data analysis. The results can be seen in the CloverGBS.html file.

SNP density
-----------

SNPdensity.R, script which displays the SNP density in all of the clover genome.

``` r
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)

snp_density <- read.csv("snp_density.csv")
snp_density$SNP.density <- as.numeric(levels(snp_density$SNP.density))[snp_density$SNP.density]

convert_ranges <- function(x) {
    values <- as.numeric(unlist(strsplit(x, "-")))
    mean(values)
}

convert_ranges <- Vectorize(convert_ranges)

binsize <- 100000

snp_density$Location <- convert_ranges(levels(snp_density$Location)[snp_density$Location])
snp_density$Coverage <- snp_density$Sites/binsize

chromosomes <- sort(as.vector(unique(snp_density[,1])))

cat("")

snp_density %>% filter(Coverage>0.01) -> snp_density

max_location <- max(snp_density$Location)

cat("Number of Sites:", sum(snp_density$Sites))
cat("Number of SNPs:", sum(snp_density$SNPs))
cat("Mean SNP density:", mean(snp_density$SNP.density))
cat("Minimum SNP density:", min(snp_density$SNP.density))
cat("Maximum SNP density:", max(snp_density$SNP.density))
cat("Variance of SNP density:", var(snp_density$SNP.density))

figures <- list()
for(chromosome in chromosomes[-1]){
    chrsub <- snp_density[snp_density$Chromosome==chromosome,]
    fig <- ggplot(chrsub, aes(x=Location, y=SNP.density)) + geom_bar(stat="identity") +
        theme_classic() + labs(title=chromosome) +
        theme(text = element_text(size=6)) +
        scale_y_continuous(limits=c(0, 0.2)) +
        scale_x_continuous(limits=c(0, max_location))
    figures[[chromosome]] <- fig
}

grid.arrange(grobs = figures, ncol=1)
```

Site Frequency Spectrum (SFS)
-----------------------------

SFS.R, script which displays the Site Frequency Spectrum (SFS).

``` r
library(ggplot2)

spectrum <- read.csv("sfs.csv")
spectrum <- spectrum[2:nrow(spectrum),1:2]

cat("Number of sites:", sum(spectrum$Count))

cat("MAF sites:", sum(spectrum$Count[1:20]))

cat("After removing MAF:", sum(spectrum$Count[-1:-20]))

spectrum$Proportion <- spectrum$Count/sum(spectrum$Count)

#spectrum <- spectrum[1:(floor(nrow(spectrum)/2)-1), 2]+spectrum[nrow(spectrum):(floor(nrow(spectrum)/2)+2), 2]

expected <- (1/1:393)/sum(1/1:393)

spectrum$Proportion[1:(ceiling(393/2)+1)] <- spectrum$Proportion[1:(ceiling(393/2)+1)] + spectrum$Proportion[394:(ceiling(393/2))]

expected <- expected[1:(ceiling(393/2)+1)]+expected[394:(ceiling(393/2))]

spectrum <- spectrum[1:(ceiling(393/2)+1), 1:3]


ggplot(spectrum, aes(x=AN, y=Proportion)) +
    geom_bar(stat="identity", fill="cyan", color="black", size=0.2) +
    labs(y="Proportion", x="Allele Frequency") + theme_classic() +
    geom_line(data=NULL, aes(x=1:198, y=expected), size=0.7, alpha=0.7)


#qplot(as.vector(rep(1:198, spectrum)), geom="histogram", binwidth=1,
#      color=I("black"), fill=I("cyan")) + theme_classic() + labs(y="Count", x="Allele Frequency")
```

Genetic Relationship Matrix (GRM)
---------------------------------

DispalyGRM.R, displays the GRM, and PCAs of the GRM, and removes the weird outliers manually.

``` r
library(ggplot2)
library(reshape2)
library(dplyr)

GRM <- as.matrix(read.csv("GRM.csv", header=FALSE))

Plants <- read.csv("White_Clover_Individuals_Library_Info.csv")

cat("Dim of GRM", dim(GRM))

seq_names <- c('1', '10', '100', '101', '102', '104', '105', '106', '107', '108', '109', '11', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '12', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '13', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '14', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '15', '150', '151', '152', '153', '154', '155', '156', '157', '158', '159', '16', '160', '161', '162', '163', '164', '165', '166', '167', '168', '169', '17', '170', '171', '172', '173', '174', '175', '176', '177', '178', '179', '18', '180', '181', '182', '183', '184', '185', '186', '187', '188', '189', '19', '190', '191', '192', '193', '194', '195', '196', '197', '198', '199', '2', '20', '200', '201', '21', '22', '23', '25', '26', '27', '28', '29', '3', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '4', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '5', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '6', '60', '61', '62', '64', '65', '66', '67', '69', '7', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '8', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '9', '90', '91', '92', '94', '95', '96', '97', '98', '99')

cat("Length of seq_names", length(seq_names))

rownames(GRM) <- seq_names
colnames(GRM) <- seq_names

GRM <- GRM[as.character(sort(as.integer(rownames(GRM)))),
           as.character(sort(as.integer(colnames(GRM))))]

collectFastaName <- Vectorize(function(x) as.integer(unlist(strsplit(as.character(x), "\\."))[1]))

plantnames <- collectFastaName(Plants$fastq.file)

indices <- vector()
for(n in rownames(GRM)){
    indices <- c(indices, which(plantnames==n))
}

all_names <- as.character(Plants$Variety[indices])
for(i in seq_along(all_names)){
    all_names[i] <- paste(unlist(strsplit(all_names[i], "\\_"))[2:3], collapse ="_")
}

rownames(GRM) <- all_names
colnames(GRM) <- all_names

GRM <- GRM[as.character(sort(rownames(GRM))),
           as.character(sort(colnames(GRM)))]

#rownames(GRM) <- 1:nrow(GRM)
#colnames(GRM) <- 1:nrow(GRM)

GRM[GRM>1] <- 1

GRMdata <- melt(GRM)

##heatmap(GRM)

ggplot(data=GRMdata, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme_classic() +
    scale_fill_gradientn(colours=c("black", "purple", "pink", "white")) +
    labs(x = "Individual", y = "Individual") + theme(axis.text.x = element_text(size=4),
                                                     axis.text.y = element_text(size=4))

## PCA

colorvector <- c("#a76d00","#4168ff","#59e069","#420097","#acd44f",
                 "#ff71f1","#00931d","#de006b","#47dada","#fe0e16",
                 "#002552","#9a9900","#70006e","#2d6800","#980016",
                 "#003003","#f6b89c","#4e0726","#ab9173","#504200")

collectGroup <- Vectorize(function(x) unlist(strsplit(as.character(x), "\\_"))[1])

PCA <- princomp(GRM, scores=TRUE)

plot(PCA)

plotdata <- data.frame(PC1=PCA$scores[,1],
                       PC2=PCA$scores[,2],
                       Groups=collectGroup(names(PCA$scores[,1])))

find_hull <- function(df) df[chull(df$PC1, df$PC2), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC1, Subdf$PC2),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC1, y=PC2, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC1, y=PC2, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 1", y="Principle component 2")


plotdata <- data.frame(PC2=PCA$scores[,2],
                       PC3=PCA$scores[,3],
                       Groups=collectGroup(names(PCA$scores[,2])))

find_hull <- function(df) df[chull(df$PC2, df$PC3), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC2, Subdf$PC3),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC2, y=PC3, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC2, y=PC3, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 2", y="Principle component 3")


plotdata <- data.frame(PC1=PCA$scores[,1],
                       PC3=PCA$scores[,3],
                       Groups=collectGroup(names(PCA$scores[,1])))

find_hull <- function(df) df[chull(df$PC1, df$PC3), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC1, Subdf$PC3),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC1, y=PC3, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC1, y=PC3, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 1", y="Principle component 3")



## Fixing mismatches

## Manually fixing labels

GRM <- GRM[-which(rownames(GRM)=="Aberconcor_07"), -which(rownames(GRM)=="Aberconcor_07")]
GRM <- GRM[-which(rownames(GRM)=="Aberconcor_03"), -which(rownames(GRM)=="Aberconcor_03")]
GRM <- GRM[-which(rownames(GRM)=="Aberpearl_10"), -which(rownames(GRM)=="Aberpearl_10")]
GRM <- GRM[-which(rownames(GRM)=="Avalon_05"), -which(rownames(GRM)=="Avalon_05")]
GRM <- GRM[-which(rownames(GRM)=="Avalon_06"), -which(rownames(GRM)=="Avalon_06")]
GRM <- GRM[-which(rownames(GRM)=="Barblanca_01"), -which(rownames(GRM)=="Barblanca_01")]
GRM <- GRM[-which(rownames(GRM)=="Borek_08"), -which(rownames(GRM)=="Borek_08")]
GRM <- GRM[-which(rownames(GRM)=="Brianna_01"), -which(rownames(GRM)=="Brianna_01")]

rownames(GRM)[rownames(GRM)=="Coolfin_04"] <- "Riesling_11"
colnames(GRM)[colnames(GRM)=="Coolfin_04"] <- "Riesling_11"

GRM <- GRM[as.character(sort(rownames(GRM))),
           as.character(sort(colnames(GRM)))]

cat("Dim of GRM: ", dim(GRM))

GRMdata <- melt(GRM)

ggplot(data=GRMdata, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme_classic() +
    scale_fill_gradientn(colours=c("black", "purple", "pink", "white")) +
    labs(x = "Individual", y = "Individual") + theme(axis.text.x = element_text(size=4),
                                                     axis.text.y = element_text(size=4))

collectGroup <- Vectorize(function(x) unlist(strsplit(as.character(x), "\\_"))[1])

PCA <- princomp(GRM, scores=TRUE)

plot(PCA)

plotdata <- data.frame(PC1=PCA$scores[,1],
                       PC2=PCA$scores[,2],
                       Groups=collectGroup(names(PCA$scores[,1])))

find_hull <- function(df) df[chull(df$PC1, df$PC2), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC1, Subdf$PC2),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC1, y=PC2, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC1, y=PC2, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 1", y="Principle component 2")


plotdata <- data.frame(PC2=PCA$scores[,2],
                       PC3=PCA$scores[,3],
                       Groups=collectGroup(names(PCA$scores[,2])))

find_hull <- function(df) df[chull(df$PC2, df$PC3), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC2, Subdf$PC3),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC2, y=PC3, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC2, y=PC3, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 2", y="Principle component 3")


plotdata <- data.frame(PC1=PCA$scores[,1],
                       PC3=PCA$scores[,3],
                       Groups=collectGroup(names(PCA$scores[,1])))

find_hull <- function(df) df[chull(df$PC1, df$PC3), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC1, Subdf$PC3),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC1, y=PC3, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC1, y=PC3, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 1", y="Principle component 3")
```

Heterozygosity
--------------

Modifcation of the SNP density, where all of the SNP densities have been divided by the 'length of species tree', or the sum of 1/1+1/2...1/(2n-1)

``` r
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)

snp_density <- read.csv("snp_density.csv")
snp_density$SNP.density <- as.numeric(levels(snp_density$SNP.density))[snp_density$SNP.density]

convert_ranges <- function(x) {
    values <- as.numeric(unlist(strsplit(x, "-")))
    mean(values)
}

convert_ranges <- Vectorize(convert_ranges)

binsize <- 10000

snp_density$Location <- convert_ranges(levels(snp_density$Location)[snp_density$Location])
snp_density$Coverage <- snp_density$Sites/binsize

chromosomes <- sort(as.vector(unique(snp_density[,1])))

snp_density %>% filter(Coverage>0.05) -> snp_density

max_location <- max(snp_density$Location)

cat("Number of Sites:", sum(snp_density$Sites))
cat("Number of SNPs:", sum(snp_density$SNPs))

snp_density$Heterozygosity <- snp_density$SNP.density/(sum(1/1:399))

cat("Mean Heterozygosity:", mean(snp_density$Heterozygosity))

figures <- list()
for(chromosome in chromosomes[-1]){
    chrsub <- snp_density[snp_density$Chromosome==chromosome,]
    fig <- ggplot(chrsub, aes(x=Location, y=Heterozygosity)) + geom_bar(stat="identity") +
        theme_classic() + labs(title=chromosome) +
        theme(text = element_text(size=6)) +
        scale_y_continuous(limits=c(0, 0.05)) +
        scale_x_continuous(limits=c(0, max_location))
    figures[[chromosome]] <- fig
}

grid.arrange(grobs = figures, ncol=1)
```

PSMC Pipeline
=============

The pipeline for Clover Full sequencing of 4 individuals used for the PSMC analysis.

gwf workflow for the PSMC

``` python
from gwf import Workflow

gwf = Workflow()

def bwa_map(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '4g',
        'cores': '8',
        'walltime': '240:00:00',
        'account': 'NChain'
    }

    spec = "./map.sh "+" ".join(infiles+outfiles)

    return options, spec

gwf.target("NCL08") << bwa_map(["NCL-08/FCHGMFMBBXX-wHAXPI040514-50_L8_1.fq.gz",
                                "NCL-08/FCHGMFMBBXX-wHAXPI040514-50_L8_2.fq.gz",
                                "NCL-08/FCHGMFNBBXX-wHAXPI040514-50_L2_1.fq.gz",
                                "NCL-08/FCHGMFNBBXX-wHAXPI040514-50_L2_2.fq.gz"],
                               ["ncl-08L8.bam", "ncl-08L2.bam"])

gwf.target("NCL09") << bwa_map(["NCL-09/FCHGMFMBBXX-wHAXPI040513-52_L8_1.fq.gz",
                                "NCL-09/FCHGMFMBBXX-wHAXPI040513-52_L8_2.fq.gz",
                                "NCL-09/FCHGMFNBBXX-wHAXPI040513-52_L2_1.fq.gz",
                                "NCL-09/FCHGMFNBBXX-wHAXPI040513-52_L2_2.fq.gz"],
                               ["ncl-09L8.bam", "ncl-09L2.bam"])

gwf.target("NCL10") << bwa_map(["NCL-10/FCHGMFMBBXX-wHAXPI040512-53_L8_1.fq.gz",
                                "NCL-10/FCHGMFMBBXX-wHAXPI040512-53_L8_2.fq.gz",
                                "NCL-10/FCHGMFNBBXX-wHAXPI040512-53_L2_1.fq.gz",
                                "NCL-10/FCHGMFNBBXX-wHAXPI040512-53_L2_2.fq.gz"],
                               ["ncl-10L8.bam", "ncl-10L2.bam"])

gwf.target("NCL12") << bwa_map(["NCL-12/FCHGMFMBBXX-wHAXPI040511-54_L8_1.fq.gz",
                                "NCL-12/FCHGMFMBBXX-wHAXPI040511-54_L8_2.fq.gz",
                                "NCL-12/FCHGMFNBBXX-wHAXPI040511-54_L2_1.fq.gz",
                                "NCL-12/FCHGMFNBBXX-wHAXPI040511-54_L2_2.fq.gz"],
                               ["ncl-12L8.bam", "ncl-12L2.bam"])


def merge_reads(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '4g',
        'walltime': '144:00:00',
        'account': 'NChain'
    }

    spec = "./merge.sh "+" ".join(infiles+outfiles)

    return options, spec

gwf.target("Merge_NCL08") << merge_reads(["ncl-08L8.bam", "ncl-08L2.bam"],
                                         ["ncl-08.u.bam"])

gwf.target("Merge_NCL09") << merge_reads(["ncl-09L8.bam", "ncl-09L2.bam"],
                                         ["ncl-09.u.bam"])

gwf.target("Merge_NCL10") << merge_reads(["ncl-10L8.bam", "ncl-10L2.bam"],
                                         ["ncl-10.u.bam"])

gwf.target("Merge_NCL12") << merge_reads(["ncl-12L8.bam", "ncl-12L2.bam"],
                                         ["ncl-12.u.bam"])

def sort_reads(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '4g',
        'walltime': '144:00:00',
        'account': 'NChain'
    }

    spec = "./sort.sh "+" ".join(infiles+outfiles)

    return options, spec

gwf.target("Sort_NCL08") << sort_reads(["ncl-08.u.bam"],
                                       ["ncl-08.bam"])

gwf.target("Sort_NCL09") << sort_reads(["ncl-09.u.bam"],
                                       ["ncl-09.bam"])

gwf.target("Sort_NCL10") << sort_reads(["ncl-10.u.bam"],
                                       ["ncl-10.bam"])

gwf.target("Sort_NCL12") << sort_reads(["ncl-12.u.bam"],
                                       ["ncl-12.bam"])


def indexbam(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '1g',
        'walltime': '144:00:00',
        'account': 'NChain'
    }

    spec = "./indexbam.sh "+" ".join(infiles+outfiles)

    return options, spec

gwf.target("Index_NCL08") << indexbam(["ncl-08.bam"],
                                      ["ncl-08.bai"])

gwf.target("Index_NCL09") << indexbam(["ncl-09.bam"],
                                      ["ncl-09.bai"])

gwf.target("Index_NCL10") << indexbam(["ncl-10.bam"],
                                      ["ncl-10.bai"])

gwf.target("Index_NCL12") << indexbam(["ncl-12.bam"],
                                      ["ncl-12.bai"])

def bamToDiploid(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '1g',
        'walltime': '72:00:00',
        'account': 'NChain'
    }

    spec = "./bam2diploid.sh "+" ".join(infiles+outfiles)

    return options, spec

gwf.target("diploid_NCL08") << bamToDiploid(["ncl-08.bam"],
                                            ["ncl-08.fq.gz"])

gwf.target("diploid_NCL09") << bamToDiploid(["ncl-09.bam"],
                                            ["ncl-09.fq.gz"])

gwf.target("diploid_NCL10") << bamToDiploid(["ncl-10.bam"],
                                            ["ncl-10.fq.gz"])

gwf.target("diploid_NCL12") << bamToDiploid(["ncl-12.bam"],
                                            ["ncl-12.fq.gz"])

def PSMC(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '1g',
        'walltime': '24:00:00',
        'account': 'NChain'
    }

    spec = "./psmc.sh "+" ".join(infiles)
    
    return options, spec

gwf.target("psmc_NCL08") << PSMC(["ncl-08.fq.gz"],
                                 ["ncl-08.fq.psmc", "ncl-08.fq.psmcfa"])

gwf.target("psmc_NCL09") << PSMC(["ncl-09.fq.gz"],
                                 ["ncl-09.fq.psmc", "ncl-09.fq.psmcfa"])

gwf.target("psmc_NCL10") << PSMC(["ncl-10.fq.gz"],
                                 ["ncl-10.fq.psmc", "ncl-10.fq.psmcfa"])

gwf.target("psmc_NCL12") << PSMC(["ncl-12.fq.gz"],
                                 ["ncl-12.fq.psmc", "ncl-12.fq.psmcfa"])


def CombinePSMC(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '1g',
        'walltime': '1:00:00',
        'account': 'NChain'
    }

    spec = "cat "+" ".join(infiles)+">> "+outfiles[0]
    
    return options, spec


gwf.target("combined") << CombinePSMC(["ncl-08.fq.psmc",
                                       "ncl-09.fq.psmc",
                                       "ncl-10.fq.psmc",
                                       "ncl-12.fq.psmc"],
                                     ["combined.psmc"])

def PlotPSMC(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '1g',
        'walltime': '24:00:00',
        'account': 'NChain'
    }

    spec = "./plotpsmc.sh "+" ".join(infiles)

    return options, spec

gwf.target("plotPSMC") << PlotPSMC(["combined.psmc"],
                                   ["combined.eps"])
```

map.sh, script for mapping reads using BWA

``` python
#!/bin/bash

source /com/extra/bwa/0.7.5a/load.sh
source /com/extra/samtools/1.3/load.sh

bwa mem -t 8 -R '@RG\tID:1\tSM:1' /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
    $1 $2 | samtools view -Sb - > $5

bwa mem -t 8 -R '@RG\tID:1\tSM:1' /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
    $3 $4 | samtools view -Sb - > $6
```

merge.sh, script for merging two bam files together

``` python
#!/bin/bash

source /com/extra/samtools/1.3/load.sh

samtools merge $3 $1 $2
```

merge.sh, script for sorting a bam file

``` python
#!/bin/bash

source /com/extra/samtools/1.3/load.sh
source /com/extra/java/8/load.sh

samtools sort -o $2 $1
```

indexbam.sh, script for indexing a bam file

``` python
#!/bin/bash

source /com/extra/samtools/1.3/load.sh

samtools index $1

mv $1.bai $2
```

bam2diploid.sh, script for coverting a bam file into a diploid fasta file.

``` r
#!/bin/bash

source /com/extra/samtools/1.4.1/load.sh
source /com/extra/bcftools/1.4.1/load.sh

samtools mpileup -C 50 -uf /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta $1 | bcftools call -c | vcfutils.pl vcf2fq -d 5 -D 100 - | gzip > ${1%.*}.fq.gz
```

psmc.sh, script for doing the psmc analysis.

``` r
#!/bin/bash

source /com/extra/psmc/2012-11-19/load.sh

fq2psmcfa $1 > ${1%.*}.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${1%.*}.psmc ${1%.*}.psmcfa
```

plotpsmc.sh, script for plotting the combined psmc results.

``` r
#!/bin/bash

source /com/extra/psmc/2012-11-19/load.sh

psmc_plot.pl -u 7e-9 -g 1 combined combined.psmc
```

Simulations
-----------

The simulations of the different hypothesis were done using a python package called msprime.

The simulation script

``` python
import msprime
from optparse import OptionParser

def help():
    print "====== Simulation through msprime with a simple growing population ====="
    print ""
    print "-o <filename.vcf>                           Output vcf file"
    print "-N int                                      Effective population size"
    print "-m float                                    Mutation rate"
    print "-r float                                    Recombination rate"
    print "-l int                                      Genome length"
    print "-g float                                    Growth rate"
    print "-s float                                    Sample size"
    print "-t int                                      Time"
    print "--print                                     Printing debug info"
    exit()

if __name__=="__main__":

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename.vcf>")
    parser.add_option('-N', type="float", nargs=1, dest="Popsize", help="Effective population size")
    parser.add_option('-m', type="float", nargs=1, dest="mutrate", help="Mutation rate")
    parser.add_option('-r', type="float", nargs=1, dest="recomrate", help="Recombination rate")
    parser.add_option('-l', type="float", nargs=1, dest="length", help="Length")
    parser.add_option('-g', type="float", nargs=1, dest="growth", help="growth rate")
    parser.add_option('-s', type="int", nargs=1, dest="samples", help="Sample size")
    parser.add_option('-t', type="float", nargs=1, dest="time", help="Time")
    parser.add_option('-d', type="int", nargs=1, dest="timeoff", help="Time of population decreasion")
    parser.add_option('-b', type="int", nargs=1, dest="bottleneck", help="Bottleneck size")
    parser.add_option('--print', action="store_true", dest="prints", help="Print debug")
    parser.add_option('--start', action="store_true", dest="start", help="Trigger Population start with bottleneck")
    parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
    options, args = parser.parse_args()

    if options.help!=None:
        help()

    if options.Popsize!=None:
        Ne = options.Popsize
    else:
        Ne = 10000
    if options.mutrate!=None:
        mu = options.mutrate
    else:
        print "No mutation rate given. This will produce no variants."
        print "use -m to add a mutation rate"
        mu = None
    if options.recomrate!=None:
        r = options.recomrate
    else:
        r = None
    if options.length!=None:
        L = options.length
    else:
        L = None
    if options.growth!=None:
        g = options.growth
    else:
        g = 0
    if options.samples!=None:
        ss = options.samples
    else:
        ss = 1
    if options.time!=None:
        t = options.time
    else:
        t = 1
    if options.timeoff!=None:
        d = options.timeoff
    else:
        d = None
    if options.bottleneck!=None:
        bottleneck_size = options.bottleneck
    else:
        bottleneck_size = None

    population_configurations = [msprime.PopulationConfiguration(sample_size=ss,
                                                                 initial_size=Ne,
                                                                 growth_rate=0)]

    if d==None and bottleneck_size==None:
        demographic_events = [msprime.PopulationParametersChange(t, growth_rate=0)]
    elif d!=None and bottleneck_size==None:
        demographic_events = [msprime.PopulationParametersChange(t-d, growth_rate=g),
                              msprime.PopulationParametersChange(t, growth_rate=0)]
    elif d==None and bottleneck_size!=None:
        if options.start!=None:
            demographic_events = [msprime.PopulationParametersChange(t, initial_size=bottleneck_size)]
        else:
            demographic_events = [msprime.PopulationParametersChange(t, initial_size=bottleneck_size),
                                  msprime.PopulationParametersChange(t+1, initial_size=Ne)]
    else:
        raise """
Warning, using -d and -b together are not implemented. For instanstaneous bottlenecks use -b.
Otherwise use -d and -g to get the bottleneck size you want."""

    if options.prints!=None:
        dd = msprime.DemographyDebugger(Ne=Ne,
                                        population_configurations=population_configurations,
                                        demographic_events=demographic_events)
        dd.print_history()

    TS = msprime.simulate(Ne=Ne, length=L, recombination_rate=r, mutation_rate=mu,
                          population_configurations=population_configurations,
                          demographic_events=demographic_events)

    if options.output!=None:
        out = options.output
    else:
        out = "output.vcf"

    with open(out, "w") as vcf_file:
        TS.write_vcf(vcf_file, 2)
```

The help message

``` r
Usage: simulate.py [options]

Options:
  -h, --help     show this help message and exit
  -o OUTPUT      <filename.vcf>
  -N POPSIZE     Effective population size
  -m MUTRATE     Mutation rate
  -r RECOMRATE   Recombination rate
  -l LENGTH      Length
  -g GROWTH      growth rate
  -s SAMPLES     Sample size
  -t TIME        Time
  -d TIMEOFF     Time of population decreasion
  -b BOTTLENECK  Bottleneck size
  --print        Print debug
  --start        Trigger Population start with bottleneck
```

Divergent genes pipeline
========================

Divergence analysis pipeline, running a given set of divergent genes and a randomly sampled set of genes

gwf workflow for the pipeline

``` python
from gwf import Workflow

gwf = Workflow()

#######################################
## Read candidate genes
#######################################

def readgenelist(filename):
    return open(filename).read().split("\n")

def cleangenelist(genelist, title):
    nlist = []
    for gene in genelist:
        if gene=="": continue
        gene = gene.split("To")[1]
        nlist.append(gene)
    return nlist

genelist = readgenelist("Divergentgenes.txt")
occi_gl = cleangenelist(genelist, "occidentale")

#############################################
## GMAP database
#############################################

def split_subgenomes(reference, subgenome1, subgenome2):
    """ """
    inputs = [reference]
    outputs = [subgenome1, subgenome2]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    directory = reference.split("/")[0]

    spec = '''
    source activate python2

    python splitreference.py {} {} {} '{}'
    '''.format(reference, subgenome1, subgenome2, ">chr9")

    return inputs, outputs, options, spec

gwf.target_from_template("RepensSplit",
                         split_subgenomes("repens/TrR.v5.fasta",
                                          "gmap/TrR.v5.To.fasta",
                                          "gmap/TrR.v5.Tp.fasta"))


#############################################
## GMAP database
#############################################

def GMAP_prepare(reference, dbname):
    """ """
    inputs = [reference]
    outputs = [dbname]
    options = {
        'cores': 1,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '04:00:00'
    }

    directory = dbname.split("/")[0]

    spec = '''
    source /com/extra/gmap/2017-08-15/load.sh

    gmap_build -d {} {} -D /home/marnit/NChain/faststorage/WHITE_CLOVER/GENE_DIVERGENCE/{}
    '''.format(dbname.split("/")[-1], reference, directory)

    return inputs, outputs, options, spec


gwf.target_from_template("RepensToDB",
                         GMAP_prepare("gmap/TrR.v5.To.fasta", "gmap/TrR.v5.To.gmap"))
gwf.target_from_template("RepensTpDB",
                         GMAP_prepare("gmap/TrR.v5.Tp.fasta", "gmap/TrR.v5.Tp.gmap"))
gwf.target_from_template("PallescensDB",
                         GMAP_prepare("pallescens/final_pallescens.fa", "gmap/pallescens.gmap"))

#############################################
## GMAP mapping on Repens and Pallescenes
#############################################

def GMAP_mapping(dbname, genes, gfffile):
    """ """
    inputs = [dbname, genes]
    outputs = [gfffile]
    options = {
        'cores': 4,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '24:00:00'
    }

    directory = dbname.split("/")[0]

    spec = '''
    source /com/extra/gmap/2017-08-15/load.sh

    gmap -d {} -D /home/marnit/NChain/faststorage/WHITE_CLOVER/GENE_DIVERGENCE/{} -t 4 -f 2 {} > {}

    '''.format(dbname.split("/")[1], directory, genes, gfffile, gfffile, gfffile)

    return inputs, outputs, options, spec


gwf.target_from_template("RepensToMapping",
                         GMAP_mapping("gmap/TrR.v5.To.gmap", "occidentale/clover.cds.fa", "gmap/TrR.To.gmap.gff"))
gwf.target_from_template("RepensTpMapping",
                         GMAP_mapping("gmap/TrR.v5.Tp.gmap", "occidentale/clover.cds.fa", "gmap/TrR.Tp.gmap.gff"))
gwf.target_from_template("PallescensMapping",
                         GMAP_mapping("gmap/pallescens.gmap", "occidentale/clover.cds.fa",
                                      "gmap/pallescens.gmap.gff"))


#############################################
## Collect all of the genes
#############################################

def collect_all_genes(inputfile, outputfile, genomefile):
    """ """
    inputs = [inputfile, genomefile]
    outputs = [outputfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '00:20:00'
    }

    spec = '''
    source /com/extra/cufflinks/2.2.1/load.sh

    gffread {} -g {} -x {}
    '''.format(inputfile, genomefile, outputfile)

    return inputs, outputs, options, spec

gwf.target_from_template("RepenToGenes",
                         collect_all_genes("gmap/TrR.To.gmap.gff", "gmap/TrR.To.gmap.fa", "gmap/TrR.v5.To.fasta"))
gwf.target_from_template("RepenTpGenes",
                         collect_all_genes("gmap/TrR.Tp.gmap.gff", "gmap/TrR.Tp.gmap.fa", "gmap/TrR.v5.Tp.fasta"))
gwf.target_from_template("PallescensGenes",
                         collect_all_genes("gmap/pallescens.gmap.gff", "gmap/pallescens.gmap.fa",
                                           "pallescens/final_pallescens.fa"))


#############################################
## Generate isolate gene files in occidentale
#############################################

def filter_gene(inputfile, outputfile, searchterm):
    """ """
    inputs = [inputfile]
    outputs = [outputfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '00:20:00'
    }

    spec = '''
    source activate python2

    python FASTafilter.py '{}' {} > {}
    '''.format(searchterm, inputfile, outputfile)

    return inputs, outputs, options, spec

for gene in occi_gl:
    gwf.target_from_template("FilterTo%s" % gene,
                             filter_gene("occidentale/clover.cds.fa",
                                        "genes/To{}.fa".format(gene),
                                        gene))
    gwf.target_from_template("FilterTp%s" % gene,
                             filter_gene("gmap/pallescens.gmap.fa",
                                        "genes/Tp{}.fa".format(gene),
                                        gene))
    gwf.target_from_template("FilterTrTo%s" % gene,
                             filter_gene("gmap/TrR.To.gmap.fa",
                                        "genes/TrTo{}.fa".format(gene),
                                        gene))
    gwf.target_from_template("FilterTrTp%s" % gene,
                             filter_gene("gmap/TrR.Tp.gmap.fa",
                                        "genes/TrTp{}.fa".format(gene),
                                        gene))

#########################################
## Multiple alignment of the genes
#########################################

def join_genes(sequences, outfile):
    inputs = sequences
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '00:20:00'
    }

    spec = '''
    source activate python2
    python joingenes.py'''
    for sequence in sequences:
        spec += " {} ".format(sequence)
    spec += "> {}".format(outfile)

    return inputs, outputs, options, spec


def multiple_align(genes, alignment):
    """ """
    inputs = [genes]
    outputs = [alignment+".best.fas"]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }

    spec = '''
    source /com/extra/mafft/7.245/load.sh

    mafft --maxiterate 1500 --globalpair --thread 2 {} > {}
    '''.format(genes, alignment)

    spec = '''
    source activate prank

    prank -d={} -o={}
    '''.format(genes, alignment)

    return inputs, outputs, options, spec

for gene in occi_gl:
    gwf.target_from_template('makeinput{}'.format(gene),
                             join_genes(['genes/To{}.fa'.format(gene),
                                         'genes/TrTo{}.fa'.format(gene),
                                         'genes/TrTp{}.fa'.format(gene),
                                         'genes/Tp{}.fa'.format(gene)],
                                        "joined_genes/{}.fasta".format(gene)))

for gene in occi_gl:
    gwf.target_from_template('alignment{}'.format(gene),
                             multiple_align("joined_genes/{}.fasta".format(gene),
                                            "aligned/{}".format(gene)))


#########################################
## Calculate dN/dS between the genes
#########################################

def calculate_dNds(alignment, dnds, pdist, sites):
    inputs = [alignment]
    outputs = [dnds, pdist]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }

    spec = '''
    source activate python2

    python dNdS.py {} {} {} {}
    '''.format(alignment, dnds, pdist, sites)

    return inputs, outputs, options, spec

for gene in occi_gl:
    gwf.target_from_template('dNdS{}'.format(gene),
                             calculate_dNds("aligned/{}.best.fas".format(gene),
                                            "summarydata/{}.dnds.csv".format(gene),
                                            "summarydata/{}.pdist.csv".format(gene),
                                            "summarydata/{}.sites.csv".format(gene)))

############################################
## Sample random genes and run full pipeline
############################################


## LOAD random genes file: random_genes.txt
## Was created by:
## grep '>' occidentale/clover.cds.fa | python createRandomList.py random_genes.txt

randomgenelist = readgenelist("random_genes.txt")

for gene in randomgenelist:
    gwf.target_from_template("RFilterTo%s" % gene,
                             filter_gene("occidentale/clover.cds.fa",
                                        "random_genes/To{}.fa".format(gene),
                                        gene))
    gwf.target_from_template("RFilterTp%s" % gene,
                             filter_gene("gmap/pallescens.gmap.fa",
                                        "random_genes/Tp{}.fa".format(gene),
                                        gene))
    gwf.target_from_template("RFilterTrTo%s" % gene,
                             filter_gene("gmap/TrR.To.gmap.fa",
                                        "random_genes/TrTo{}.fa".format(gene),
                                        gene))
    gwf.target_from_template("RFilterTrTp%s" % gene,
                             filter_gene("gmap/TrR.Tp.gmap.fa",
                                        "random_genes/TrTp{}.fa".format(gene),
                                        gene))

for gene in randomgenelist:
    gwf.target_from_template('Rmakeinput{}'.format(gene),
                             join_genes(['random_genes/To{}.fa'.format(gene),
                                         'random_genes/TrTo{}.fa'.format(gene),
                                         'random_genes/TrTp{}.fa'.format(gene),
                                         'random_genes/Tp{}.fa'.format(gene)],
                                        "joined_random_genes/{}.fasta".format(gene)))

for gene in randomgenelist:
    gwf.target_from_template('Ralignment{}'.format(gene),
                             multiple_align("joined_random_genes/{}.fasta".format(gene),
                                            "aligned_random/{}".format(gene)))


for gene in randomgenelist:
    gwf.target_from_template('RdNdS{}'.format(gene),
                             calculate_dNds("aligned_random/{}.best.fas".format(gene),
                                            "summarydata_random/{}.dnds.csv".format(gene),
                                            "summarydata_random/{}.pdist.csv".format(gene),
                                            "summarydata_random/{}.sites.csv".format(gene)))

def extract_from_vcf(gene, gff_file, vcf_file, output, outname):
    inputs = [gff_file, vcf_file]
    outputs = [output]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }

    spec = '''
    source activate python2
    source /com/extra/bedops/2.4.26/load.sh
    source /com/extra/tabix/0.2.6/load.sh

    cd vcf_data

    grep '{gene}' ../{gff} | grep 'CDS' | gff2bed | python ../tabixsearch.py ../{vcf} | python ../vcf2fasta.py {name} > ../{out}
    '''.format(gene=gene, gff=gff_file, vcf=vcf_file, out=output, name=outname)

    return inputs, outputs, options, spec

for gene in occi_gl:
    gwf.target_from_template('NCL08TrTovcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.To.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-08.g.vcf.gz",
                                              "vcf_data/ncl-08-TrTo"+gene+".fasta", "ncl-08-TrTo"+gene))
    gwf.target_from_template('NCL08TrTpvcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.Tp.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-08.g.vcf.gz",
                                              "vcf_data/ncl-08-TrTp"+gene+".fasta", "ncl-08-TrTp"+gene))

    gwf.target_from_template('NCL09TrTovcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.To.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-09.g.vcf.gz",
                                              "vcf_data/ncl-09-TrTo"+gene+".fasta", "ncl-09-TrTo"+gene))
    gwf.target_from_template('NCL09TrTpvcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.Tp.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-09.g.vcf.gz",
                                              "vcf_data/ncl-09-TrTp"+gene+".fasta", "ncl-09-TrTp"+gene))

    gwf.target_from_template('NCL10TrTovcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.To.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-10.g.vcf.gz",
                                              "vcf_data/ncl-10-TrTo"+gene+".fasta", "ncl-10-TrTo"+gene))
    gwf.target_from_template('NCL10TrTpvcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.Tp.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-10.g.vcf.gz",
                                              "vcf_data/ncl-10-TrTp"+gene+".fasta", "ncl-10-TrTp"+gene))

    gwf.target_from_template('NCL12TrTovcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.To.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-12.g.vcf.gz",
                                              "vcf_data/ncl-12-TrTo"+gene+".fasta", "ncl-12-TrTo"+gene))

    gwf.target_from_template('NCL12TrTpvcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.Tp.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-12.g.vcf.gz",
                                              "vcf_data/ncl-12-TrTp"+gene+".fasta", "ncl-12-TrTp"+gene))

for gene in occi_gl:
    gwf.target_from_template('join{}'.format(gene),
                             join_genes(["vcf_data/ncl-08-TrTo"+gene+".fasta", "vcf_data/ncl-09-TrTo"+gene+".fasta",
                                         "vcf_data/ncl-10-TrTo"+gene+".fasta", "vcf_data/ncl-12-TrTo"+gene+".fasta",
                                         "vcf_data/ncl-08-TrTp"+gene+".fasta", "vcf_data/ncl-09-TrTp"+gene+".fasta",
                                         "vcf_data/ncl-10-TrTp"+gene+".fasta", "vcf_data/ncl-12-TrTp"+gene+".fasta"],
                                        "joined_vcf_data/{}.vcf.fasta".format(gene)))

    gwf.target_from_template('vcf_alignment{}'.format(gene),
                             multiple_align("joined_vcf_data/{}.vcf.fasta".format(gene),
                                            "vcf_aligned/{}".format(gene)))
for gene in occi_gl:
    gwf.target_from_template('vcf_dNdS{}'.format(gene),
                             calculate_dNds("vcf_aligned/{}.best.fas".format(gene),
                                            "vcf_summarydata/{}.dnds.csv".format(gene),
                                            "vcf_summarydata/{}.pdist.csv".format(gene),
                                            "vcf_summarydata/{}.sites.csv".format(gene)))

def summary_files(inputfiles, output):
    inputs = inputfiles
    outputs = [output]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }

    inputstring = ",".join(inputfiles)

    spec = '''
    source activate python2

    python Joinsummary.py "{}" {}
    '''.format(inputstring, output)

    return inputs, outputs, options, spec

blacklist = ["1161g41190.1", "96g69060.1", "171g53460.1", "971g00020.1"]

dNdSfiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    dNdSfiles.append("summarydata/{}.dnds.csv".format(gene))

pdistfiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    pdistfiles.append("summarydata/{}.pdist.csv".format(gene))

sitefiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    sitefiles.append("summarydata/{}.sites.csv".format(gene))

gwf.target_from_template("dNdSsummarybetween",
                         summary_files(dNdSfiles, "summaryfiles/dNdSsummaryBetween.csv"))
gwf.target_from_template("pdistsummarybetween",
                         summary_files(pdistfiles, "summaryfiles/pdistsummaryBetween.csv"))
#gwf.target_from_template("sitesummarybetween",
#                         summary_files(sitefiles, "summaryfiles/sitesummaryBetween.csv"))


blacklist = ["685g22150.1", "2084g01010.1", "6044g00010.1", "2155g00010.1", "85g75560.1"]

RdNdSfiles = []
for gene in randomgenelist:
    if gene in blacklist: continue
    RdNdSfiles.append("summarydata_random/{}.dnds.csv".format(gene))

Rpdistfiles = []
for gene in randomgenelist:
    if gene in blacklist: continue
    Rpdistfiles.append("summarydata_random/{}.pdist.csv".format(gene))

Rsitefiles = []
for gene in randomgenelist:
    if gene in blacklist: continue
    Rsitefiles.append("summarydata_random/{}.sites.csv".format(gene))

gwf.target_from_template("RdNdSsummarybetween",
                         summary_files(RdNdSfiles, "summaryfiles/random_dNdSsummaryBetween.csv"))
gwf.target_from_template("Rpdistsummarybetween",
                         summary_files(Rpdistfiles, "summaryfiles/random_pdistsummaryBetween.csv"))
#gwf.target_from_template("Rsitesummarybetween",
#                         summary_files(Rsitefiles, "summaryfiles/random_sitesummaryBetween.csv"))

blacklist = ["96g69060.1", "1460g01050.1", "1147g10020.2", "201g20030.1", "2483g10120.1",
             "3506g10030.1"]

vcfdNdSfiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    vcfdNdSfiles.append("vcf_summarydata/{}.dnds.csv".format(gene))

vcfpdistfiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    vcfpdistfiles.append("vcf_summarydata/{}.pdist.csv".format(gene))

vcfsitefiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    vcfsitefiles.append("vcf_summarydata/{}.sites.csv".format(gene))

gwf.target_from_template("vcfdNdSsummary",
                         summary_files(vcfdNdSfiles, "summaryfiles/vcf_dNdSsummary.csv"))
gwf.target_from_template("vcfpdistsummary",
                         summary_files(vcfpdistfiles, "summaryfiles/vcf_pdistsummary.csv"))
#gwf.target_from_template("vcfsitesummary",
#                         summary_files(vcfsitefiles, "summaryfiles/vcf_sitesummary.csv"))



for gene in randomgenelist:
    gwf.target_from_template('RNCL08TrTovcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.To.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-08.g.vcf.gz",
                                              "vcf_data_random/ncl-08-TrTo"+gene+".fasta", "ncl-08-TrTo"+gene))
    gwf.target_from_template('RNCL08TrTpvcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.Tp.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-08.g.vcf.gz",
                                              "vcf_data_random/ncl-08-TrTp"+gene+".fasta", "ncl-08-TrTp"+gene))

    gwf.target_from_template('RNCL09TrTovcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.To.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-09.g.vcf.gz",
                                              "vcf_data_random/ncl-09-TrTo"+gene+".fasta", "ncl-09-TrTo"+gene))
    gwf.target_from_template('RNCL09TrTpvcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.Tp.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-09.g.vcf.gz",
                                              "vcf_data_random/ncl-09-TrTp"+gene+".fasta", "ncl-09-TrTp"+gene))

    gwf.target_from_template('RNCL10TrTovcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.To.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-10.g.vcf.gz",
                                              "vcf_data_random/ncl-10-TrTo"+gene+".fasta", "ncl-10-TrTo"+gene))
    gwf.target_from_template('RNCL10TrTpvcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.Tp.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-10.g.vcf.gz",
                                              "vcf_data_random/ncl-10-TrTp"+gene+".fasta", "ncl-10-TrTp"+gene))

    gwf.target_from_template('RNCL12TrTovcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.To.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-12.g.vcf.gz",
                                              "vcf_data_random/ncl-12-TrTo"+gene+".fasta", "ncl-12-TrTo"+gene))

    gwf.target_from_template('RNCL12TrTpvcf{}'.format(gene),
                             extract_from_vcf(gene, "gmap/TrR.Tp.gmap.gff",
                                              "../CLOVER_RESEQUENCING/ncl-12.g.vcf.gz",
                                              "vcf_data_random/ncl-12-TrTp"+gene+".fasta", "ncl-12-TrTp"+gene))


for gene in randomgenelist:
    gwf.target_from_template('Rjoin{}'.format(gene),
                             join_genes(["vcf_data_random/ncl-08-TrTo"+gene+".fasta", "vcf_data_random/ncl-09-TrTo"+gene+".fasta",
                                         "vcf_data_random/ncl-10-TrTo"+gene+".fasta", "vcf_data_random/ncl-12-TrTo"+gene+".fasta",
                                         "vcf_data_random/ncl-08-TrTp"+gene+".fasta", "vcf_data_random/ncl-09-TrTp"+gene+".fasta",
                                         "vcf_data_random/ncl-10-TrTp"+gene+".fasta", "vcf_data_random/ncl-12-TrTp"+gene+".fasta"],
                                        "joined_vcf_data_random/{}.vcf.fasta".format(gene)))

    gwf.target_from_template('Rvcf_alignment{}'.format(gene),
                             multiple_align("joined_vcf_data_random/{}.vcf.fasta".format(gene),
                                            "vcf_aligned_random/{}".format(gene)))
for gene in randomgenelist:
    gwf.target_from_template('vcf_dNdS{}'.format(gene),
                             calculate_dNds("vcf_aligned_random/{}.best.fas".format(gene),
                                            "vcf_summarydata_random/{}.dnds.csv".format(gene),
                                            "vcf_summarydata_random/{}.pdist.csv".format(gene),
                                            "vcf_summarydata_random/{}.sites.csv".format(gene)))

blacklist = ["2155g00010.1", "1460g01050.1", "685g22150.1", "6044g00010.1", "473g62370.1",
             "2518g20040.1"]

RvcfdNdSfiles = []
for gene in randomgenelist:
    if gene in blacklist: continue
    RvcfdNdSfiles.append("vcf_summarydata_random/{}.dnds.csv".format(gene))

Rvcfpdistfiles = []
for gene in randomgenelist:
    if gene in blacklist: continue
    Rvcfpdistfiles.append("vcf_summarydata_random/{}.pdist.csv".format(gene))

Rvcfsitefiles = []
for gene in randomgenelist:
    if gene in blacklist: continue
    Rvcfsitefiles.append("vcf_summarydata_random/{}.sites.csv".format(gene))

gwf.target_from_template("RvcfdNdSsummary",
                         summary_files(RvcfdNdSfiles, "summaryfiles/vcf_random_dNdSsummary.csv"))
gwf.target_from_template("Rvcfpdistsummary",
                         summary_files(Rvcfpdistfiles, "summaryfiles/vcf_random_pdistsummary.csv"))
#gwf.target_from_template("Rvcfsitesummary",
#                         summary_files(Rvcfsitefiles, "summaryfiles/vcf_random_sitesummary.csv"))
```

Splitting the white clover reference into the two subgenomes

``` python
from sys import argv

def make_new_reference_files(filename, sub1, sub2, divider=">chr9"):
    genomes = open(filename).read().split(divider)
    f = open(sub1, "w")
    f.write(genomes[0])
    f.close()
    f = open(sub2, "w")
    f.write(">chr9"+genomes[1])
    f.close()

if __name__=="__main__":
    make_new_reference_files(argv[1], argv[2], argv[3], argv[4])
```

FASTafilter
-----------

Fasta fitler script, which finds a given sequence using the sequence name. First grep, then jump to the file location to readout the sequence.

``` python
import subprocess 
from sys import argv, exit

def findlocation(filename, search_term):
    out = subprocess.check_output("grep -b '{}' {}".format(search_term, filename), shell=True)
    return int(out.split(":")[0])

def filtersequence(filename, location):
    f = open(filename)
    f.seek(location)
    filebegun = False
    for line in f:
        if filebegun==False:
            if ">" in line:
                filebegun = True
                print line,
        else:
            if ">" in line: exit()
            print line,
            
if __name__=="__main__":
    search_term, filename = argv[1:3]
    location = findlocation(filename, search_term)
    filtersequence(filename, location)
```

Joins a list of given genes, and cleans the sequence names.

``` python
from sys import argv

def write_out_fasta(genes):
    for gene in genes:
        f = open(gene).read().split("\n", 1)[1]
        print ">"+gene.split(".")[0].split("/")[-1]
        print f,

if __name__=="__main__":
    genes = argv[1:]
    write_out_fasta(genes)
```

dNdS script
-----------

Calculates dNdS, p-distance and gives raw numbers back. It produces 3 different summary output files. dNdS, pdist and sites.

``` python
from sys import argv
from math import log

codon_map = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
             'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
             'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*',
             'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
             'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H',
             'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
             'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
             'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
             'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S',
             'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V',
             'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
             'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
             'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

def readfasta(filename):
    sequences = open(filename).read().split(">")[1:]
    sequences = [seq.split("\n", 1) for seq in sequences]
    return {seq[0]:seq[1].replace("\n", "") for seq in sequences}


def all_paths(wd, codon, gcodon):
    indices = filter(lambda x: x>-1, [w*(i+1)-1 for w, i in zip(wd, range(3))])
    pathways = []
    if len(indices)>1:
        for i in indices:
            ncodon = codon[:i]+gcodon[i]+codon[i+1:]
            b1, b2 = codon_map[codon], codon_map[ncodon]
            if b1=="*" or b2=="*":
                pathways.append(None)
                continue
            cpath = [0,0]
            if b1==b2:
                cpath[0] += 1
            else:
                cpath[1] += 1
            new_cases = wd[:i]+[False]+wd[i+1:]
            internal_paths = all_paths(new_cases, ncodon, gcodon)
            if internal_paths==None:
                pathways.append(None)
                continue
            if type(internal_paths[0])==int:
                internal_paths = [internal_paths]
            for path, i in zip(internal_paths, range(len(internal_paths))):
                if path==None:
                    pathways.append(None)
                    continue
                internal_paths[i][0] += cpath[0]
                internal_paths[i][1] += cpath[1]
                pathways.append(internal_paths[i])
    else:
        i = indices[0]
        ncodon = codon[:i]+gcodon[i]+codon[i+1:]
        b1, b2 = codon_map[codon], codon_map[ncodon]
        if b1=="*" or b2=="*":
            return None
        if b1==b2:
            return [1, 0]
        else:
            return [0, 1]
    return pathways

def count_pathways(codon1, codon2):
    s = 0
    n = 0
    which_different = []
    for b1, b2 in zip(codon1, codon2):
        if b1==b2:
            which_different.append(False)
        else:
            which_different.append(True)
    how_many_different = sum(which_different)
    pathways = all_paths(which_different, codon1, codon2)
    if type(pathways[0])==int:
        pathways = [pathways]
    c = 0
    for pathway in pathways:
        if pathway == None: continue
        s += pathway[0]
        n += pathway[1]
        c += 1
    return s/float(c), n/float(c)

def jukes_cantor(D):
    return -3/float(4)*log(1-4/float(3)*D)

def calculate_dn_ds(seq1, seq2):
    seq_n = len(seq1)
    Sd, Nd = 0, 0
    S, N = 0, 0
    i = 1
    for i in range(0, seq_n, 3):
        codon1, codon2 = seq1[i:i+3], seq2[i:i+3]
        base1, base2 = codon_map.get(codon1, "-"), codon_map.get(codon2, "-")
        #print base1, base2
        if base1=="-" or base2=="-": continue
        if codon1!=codon2:
            S += 1.5
            N += 1.5
            #print codon1, codon2
            sdnd = count_pathways(codon1, codon2)
            #print sdnd
            Sd += sdnd[0]
            Nd += sdnd[1]
    if S==0 and N==0:
        return 1
    #print "S & N", S, N
    #print "Sd & Nd", Sd, Nd
    pS = Sd/float(S)
    pN = Nd/float(N)
    #print "pS & pN", pS, pN
    dS = jukes_cantor(pS)
    dN = jukes_cantor(pN)
    #print "dS & dN", dS, dN
    if dS==0:
        return 1
    return dN/dS

def pdistance(seq1, seq2):
    n = 0
    d = 0
    for base1, base2 in zip(seq1, seq2):
        if base1=="-" or base2=="-": continue
        if base1!=base2:
            d += 1
        n += 1
    return d/float(n), (d, n)

def cleanKeyNames(keys):
    new_keys = []
    for key in keys:
        if "T"==key[2] and "T"==key[0]:
            new_keys.append(key[:4])
        elif "T"==key[0]:
            new_keys.append(key[:2])
        else:
            new_keys.append(key[:11])
    return new_keys

def prettyprint(matrix):
    keys = sorted(matrix.keys())
    for key in [""]+keys:
        print "{:5}".format(key),
    print "\n",
    for key in keys:
        print "{:4}".format(key),
        for key2 in keys:
            print "%.3f" % matrix[key][key2],
        print

def write_out_matrix(filename, matrix):
    keys = sorted(matrix.keys())
    f = open(filename, "w")
    f.write("ID")
    for key in keys:
        f.write(", {}".format(key))
    f.write("\n")
    for key in keys:
        f.write("{}".format(key))
        for key2 in keys:
            if type(matrix[key][key2])==tuple:
                f.write(", ({};{})".format(matrix[key][key2][0], matrix[key][key2][1]))
            else:
                f.write(", {}".format(matrix[key][key2]))
        f.write("\n")
    f.close()

if __name__=="__main__":
    sequences = readfasta(argv[1])
    done = []

    samples = cleanKeyNames(sequences.keys())
    for key, sample in zip(sequences.keys(), samples):
        sequences[sample] = sequences[key]
        del sequences[key]

    dNdS = {}
    pdist = {}
    sites = {}
    for seq1 in sequences:
        if seq1 not in dNdS:
            dNdS[seq1] = {}
            pdist[seq1] = {}
            sites[seq1] = {}
        for seq2 in sequences:
            dNdS[seq1][seq2] = calculate_dn_ds(sequences[seq1], sequences[seq2])
            pdist[seq1][seq2], sites[seq1][seq2] = pdistance(sequences[seq1], sequences[seq2])

    prettyprint(dNdS)
    prettyprint(pdist)

    write_out_matrix(argv[2], dNdS)
    write_out_matrix(argv[3], pdist)
    write_out_matrix(argv[4], sites)
```

Tabixsearch
-----------

Reads in a bed file, and uses tabix to extract a region from vcf.gz file.

``` python
import sys
import os
from stat import S_ISFIFO
import subprocess

def tabix_vcf(beddata, vcffile):
    call = "tabix "+vcffile
    for line in beddata:
        if line=='': continue
        line = line.split("\t")
        orientation = line[5]
        call += " {chrom}:{start}-{end}".format(chrom=line[0], start=line[1], end=line[2])
    subprocess.call(call, shell=True)
    print(orientation)
    
if __name__=="__main__":
    if S_ISFIFO(os.fstat(0).st_mode):
        bedfile = sys.stdin.readlines()
        vcffile = sys.argv[1]
        tabix_vcf(bedfile, vcffile)
    else:
        bedfile = open(sys.argv[1]).read().split("\n")
        vcffile = sys.argv[2]
        tabix_vcf(bedfile, vcffile)
        
```

vcf2fasta.py
------------

Read in vcf segment, and converts it into a diploid fasta format.

``` python
import sys
import os
from stat import S_ISFIFO

hetero_codes = {'A': {'G': 'R', 'T': 'W', 'C': 'M'},
                'G': {'A': 'R', 'T': 'K', 'C': 'S'},
                'T': {'A': 'W', 'G': 'K', 'C': 'Y'},
                'C': {'A': 'M', 'G': 'S', 'T': 'Y'}} 

def reverse_complement(sequence):
    comple = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y',
              'Y': 'R', 'S': 'W', 'W': 'S', 'K': 'M', 'M': 'K'}
    return "".join([comple.get(base, 'N') for base in sequence[::-1]])

def convert_vcf(vcffile, name):
    s = ""
    print ">"+name
    orient = ""
    for line in vcffile:
        if line=="-" or line=="+":
            orient = line.replace("\n", "")
        line = line.split("\t")
        if len(line)<9: continue
        chrom, pos, _, ref, alt = line[:5]
        genotype = line[9].split(":")[0]
        if alt==".":
            s += ref
        else:
            if len(alt)>1:
                alt1, alt2 = alt.split(",")
                s += hetero_codes[alt1][alt2]
                continue
            if genotype=="1/1":
                s += alt
            else:
                s += hetero_codes[ref][alt]
    if orient=="-": s = reverse_complement(s)
    for i in range(len(s)/60+1):
        print s[i*60:(i+1)*60]


if __name__=="__main__":
    if S_ISFIFO(os.fstat(0).st_mode):
        vcffile = sys.stdin.readlines()
        name = sys.argv[1]
        convert_vcf(vcffile, name)
    else:
        vcffile = open(sys.argv[1]).read().split("\n")
        name = sys.argv[2]
        convert_vcf(vcffile, name)
```

Joins summary files returned from dNdS into one summary file.

``` python
from sys import argv

def read_csv_files(filenames):
    full_data = {}
    full_data["ID"] = []
    for filename in filenames:
        f = open(filename)
        full_data["ID"].append(filename.split("/")[-1])
        row = {}
        for line in f:
            line = line.replace("\n", "")
            line = line.split(", ")
            if len(line)<2: continue
            if line[0]=="ID":
                header = line[1:]
                continue
            seq = line[0]
            for data, seq2 in zip(line[1:], header):
                if seq==seq2: continue
                combname = "|".join(sorted([seq,seq2]))
                if combname not in row:
                    row[combname] = data
        for item in row:
            if item not in full_data:
                full_data[item] = []
            full_data[item].append(row[item])
    return full_data

def write_out_data(filename, full_data):
    f = open(filename, "w")
    sequence_names = full_data["ID"]
    del full_data["ID"]
    keys = sorted(full_data.keys())
    nrows = len(sequence_names)
    f.write("ID")
    for key in keys:
        f.write(","+key)
    f.write("\n")
    for i in range(nrows):
        f.write(sequence_names[i])
        for key in keys:
            f.write(","+full_data[key][i])
        f.write("\n")
    f.close()

if __name__=="__main__":
    filenames = argv[1].split(",")
    combined_data = read_csv_files(filenames)
    write_out_data(argv[2], combined_data)
```
