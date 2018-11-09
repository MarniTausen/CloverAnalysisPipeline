Clover genotype analysis supplementary
================
09/11/2018 - 09:51:11

-   [Introduction](#introduction)
    -   [Unifiedgenotyper](#unifiedgenotyper)
    -   [Selectvariants](#selectvariants)
    -   [SNP density](#snp-density)
    -   [Site Frequency Spectrum](#site-frequency-spectrum)
    -   [Parse VCF](#parse-vcf)
        -   [LD analysis](#ld-analysis)
    -   [Genetic Relationship Matrix](#genetic-relationship-matrix)
-   [GBS data analysis](#gbs-data-analysis)
    -   [Site Frequency Spectrum (SFS)](#site-frequency-spectrum-sfs)
    -   [Genetic Relationship Matrix (GRM)](#genetic-relationship-matrix-grm)
    -   [Heterozygosity](#heterozygosity)
-   [PSMC Pipeline](#psmc-pipeline)
    -   [Simulations](#simulations)
-   [Divergent genes pipeline](#divergent-genes-pipeline)
    -   [FASTafilter](#fastafilter)
    -   [dNdS script](#dnds-script)
-   [Gene annotation pipeline](#gene-annotation-pipeline)

Introduction
============

This document has all of the scripts used for the Clover genotype analysis. The results are included in a separate document. Each step of the process has its own section and code in the correct order it was run. \#GBS genotyping The GBS genotyping of the 200 clover individuals. The main workflow script is a gwf pipeline (<http://gwf.readthedocs.io/en/latest/>), which keeps track of all job dependencies and submits the scripts in the correct order.

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

``` bash
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

``` bash
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

``` bash
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

``` bash
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

``` bash
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

``` bash
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

``` bash
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

``` r
====== SNP density analysis of a gVCF file =====
Reads in a gvcf file, strips it down and counts
the occurrences of SNPS versus the number sites sequenced
-i <filename.vcf>                 The input vcf file
-o <filename.csv>                 The output csv file
-b binsize                        The binsize to be calculated from
-d read depth                     Filter sites less that the minimum read depth
```

### Site Frequency Spectrum

``` bash
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

``` bash
#!/bin/bash

source /com/extra/vcftools/0.1.14/load.sh

/com/extra/vcftools/0.1.14/bin/vcftools --vcf Clover_SP_dbSNP_V1.1.MAF.vcf --geno-r2 --ld-window-bp 10000 --out LD_10k_window
```

./LDR.sh Script loading the Rscript which produces the figures.

``` r
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

``` bash
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

This section contains all of the R scripts which visualize the results from the GBS data analysis. The results can be seen in the CloverGBS.html file. \#\#SNP density SNPdensity.R, script which displays the SNP density in all of the clover genome.

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

The pipeline for Clover Full sequencing of 4 individuals used for the PSMC analysis. gwf workflow for the PSMC

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

``` bash
#!/bin/bash

source /com/extra/bwa/0.7.5a/load.sh
source /com/extra/samtools/1.3/load.sh

bwa mem -t 8 -R '@RG\tID:1\tSM:1' /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
    $1 $2 | samtools view -Sb - > $5

bwa mem -t 8 -R '@RG\tID:1\tSM:1' /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta \
    $3 $4 | samtools view -Sb - > $6
```

merge.sh, script for merging two bam files together

``` bash
#!/bin/bash

source /com/extra/samtools/1.3/load.sh

samtools merge $3 $1 $2
```

merge.sh, script for sorting a bam file

``` bash
#!/bin/bash

source /com/extra/samtools/1.3/load.sh
source /com/extra/java/8/load.sh

samtools sort -o $2 $1
```

indexbam.sh, script for indexing a bam file

``` bash
#!/bin/bash

source /com/extra/samtools/1.3/load.sh

samtools index $1

mv $1.bai $2
```

bam2diploid.sh, script for coverting a bam file into a diploid fasta file.

``` bash
#!/bin/bash

source /com/extra/samtools/1.4.1/load.sh
source /com/extra/bcftools/1.4.1/load.sh

samtools mpileup -C 50 -uf /home/marnit/NChain/faststorage/WHITE_CLOVER/WCL_IND/Reference/TrR.v5.fasta $1 | bcftools call -c | vcfutils.pl vcf2fq -d 5 -D 100 - | gzip > ${1%.*}.fq.gz
```

psmc.sh, script for doing the psmc analysis.

``` bash
#!/bin/bash

source /com/extra/psmc/2012-11-19/load.sh

fq2psmcfa $1 > ${1%.*}.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${1%.*}.psmc ${1%.*}.psmcfa
```

plotpsmc.sh, script for plotting the combined psmc results.

``` bash
#!/bin/bash

source /com/extra/psmc/2012-11-19/load.sh

psmc_plot.pl -u 7e-9 -g 1 combined combined.psmc
```

Simulations
-----------

The simulations of the different hypothesis were done using a python package called msprime. The simulation script

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

Divergence analysis pipeline, running a given set of divergent genes and a randomly sampled set of genes gwf workflow for the pipeline

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
#occi_gl = cleangenelist(genelist, "occidentale")
occi_gl = genelist
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
                                          "repens/TrR.v5.To.fasta",
                                          "repens/TrR.v5.Tp.fasta"))
#############################################
## Making blast databases?
#############################################
#############################################
## genblast
#############################################
def translateSequences(infile, outfile):
    inputs = [infile]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '04:00:00'
    }
    spec = '''
    source activate python2
    python translatefasta.py {} {}
    '''.format(infile, outfile)
    return inputs, outputs, options, spec
def genblast(query_genes, database, outfile, basename=""):
    """ """
    inputs = [query_genes, database]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '04:00:00'
    }
    if basename=="":
        basename = outfile.split("/")[-1].split(".")[0]
    spec = '''
    source activate python2
    cd genBlast
    ./genblast -p genblastg -q ../{query} -t ../{db} -c 1 -r 1 -gff -cdna -o {out}
    python ../cleanfasta.py {out}_1.1c_2.3_s1_0_16_1.DNA > {out}.cleaned.DNA
    python ../FASTafilter.py '\-R1\-' {out}.cleaned.DNA > ../{output}
    '''.format(query=query_genes, db=database, out=basename, output=outfile)
    return inputs, outputs, options, spec
gwf.target_from_template("RepensTogenBlast",
                         genblast("occidentale/divergentgenes.fasta", "repens/TrR.v5.To.fasta", "genblastresults/Togenes.fasta"))
gwf.target_from_template("PalgenBlast",
                         genblast("occidentale/divergentgenes.fasta", "pallescens/final_pallescens.fa", "genblastresults/Palgenes.fasta"))
gwf.target_from_template("Paltranslate",
                         translateSequences("genblastresults/Palgenes.fasta", "genblastresults/Palgenes.prot.fasta"))
gwf.target_from_template("RepensTpgenBlast",
                         genblast("genblastresults/Palgenes.prot.fasta", "repens/TrR.v5.Tp.fasta", "genblastresults/Tpgenes.fasta"))
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
                                         "_"+gene))
    gwf.target_from_template("FilterTp%s" % gene,
                             filter_gene("genblastresults/Palgenes.fasta",
                                        "genes/Tp{}.fa".format(gene),
                                         "_"+gene))
    gwf.target_from_template("FilterTrTo%s" % gene,
                             filter_gene("genblastresults/Togenes.fasta",
                                        "genes/TrTo{}.fa".format(gene),
                                         "_"+gene))
    gwf.target_from_template("FilterTrTp%s" % gene,
                             filter_gene("genblastresults/Tpgenes.fasta",
                                        "genes/TrTp{}.fa".format(gene),
                                         "_"+gene))
#########################################
## Multiple alignment of the genes
#########################################
def join_genes(sequences, outfile, reference_gene):
    inputs = sequences
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '00:20:00'
    }
    spec = '''
    source activate biopython2
    python FASTafilter.py '{}' {} | python joingenes.py'''.format(reference_gene, "occidentale/clover.cds.fa")
    for sequence in sequences:
        spec += " {} ".format(sequence)
    spec += "> {}".format(outfile)
    return inputs, outputs, options, spec
def join_genes_and_translate(sequences, outfile, reference_gene):
    inputs = sequences
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '00:20:00'
    }
    spec = '''
    source activate biopython2
    python FASTafilter.py '{}' {} | python joingenesandtranslate.py'''.format(reference_gene, "occidentale/clover.cds.fa")
    for sequence in sequences:
        spec += " {} ".format(sequence)
    spec += "> {}".format(outfile)
    return inputs, outputs, options, spec
def multiple_align(genes, alignment, settings):
    """ """
    inputs = [genes]
    outputs = [alignment]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }
    spec = '''
    source activate prank
    prank -d={} -o={} {}
    '''.format(genes, ".".join(alignment.split(".")[0:2]), settings)
    return inputs, outputs, options, spec
for gene in occi_gl:
    gwf.target_from_template('makeinput{}'.format(gene),
                             join_genes(['genes/To{}.fa'.format(gene),
                                         'genes/TrTo{}.fa'.format(gene),
                                         'genes/TrTp{}.fa'.format(gene),
                                         'genes/Tp{}.fa'.format(gene)],
                                        "joined_genes/{}.fasta".format(gene), gene))
for gene in occi_gl:
    gwf.target_from_template('alignment{}'.format(gene),
                             multiple_align("joined_genes/{}.fasta".format(gene),
                                            "aligned/{}.best.fas".format(gene), "-codon +F"))
for gene in occi_gl:
    gwf.target_from_template('pepmakeinput{}'.format(gene),
                             join_genes_and_translate(['genes/To{}.fa'.format(gene),
                                         'genes/TrTo{}.fa'.format(gene),
                                         'genes/TrTp{}.fa'.format(gene),
                                         'genes/Tp{}.fa'.format(gene)],
                                        "joined_proteins/{}.fasta".format(gene), gene))
for gene in occi_gl:
    gwf.target_from_template('pepalignment{}'.format(gene),
                             multiple_align("joined_proteins/{}.fasta".format(gene),
                                            "aligned_proteins/{}.best.fas".format(gene), "-protein"))
#########################################
## Calculate dN/dS between the genes
#########################################
def calculate_dNds(alignment, dnds, pdist, sites, dnds_data=""):
    inputs = [alignment]
    outputs = [dnds, pdist, sites, dnds_data]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }
    spec = '''
    source activate python2
    python dNdS.py {} {} {} {} {}
    '''.format(alignment, dnds, pdist, sites, dnds_data)
    return inputs, outputs, options, spec
for gene in occi_gl:
   gwf.target_from_template('dNdS{}'.format(gene),
                            calculate_dNds("aligned/{}.best.fas".format(gene),
                                           "summarydata/{}.dnds.csv".format(gene),
                                           "summarydata/{}.pdist.csv".format(gene),
                                           "summarydata/{}.sites.csv".format(gene),
                                           "summarydata/{}.dnds.info.csv".format(gene)))
def comparePeptides(alignment, pdist, sites):
    inputs = [alignment]
    outputs = [pdist, sites]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }
    spec = '''
    source activate python2
    python comparePeptides.py {} {} {}
    '''.format(alignment, pdist, sites)
    return inputs, outputs, options, spec
for gene in occi_gl:
   gwf.target_from_template('cmpPep{}'.format(gene),
                            comparePeptides("aligned_proteins/{}.best.fas".format(gene),
                                           "summarydata_protein/{}.pdist.csv".format(gene),
                                           "summarydata_protein/{}.sites.csv".format(gene)))
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
blacklist = ["840g41210.1", "0g01140.1", "4488g00060.1", "1089g30080.1", ""]
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
dNdSinfofiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    dNdSinfofiles.append("summarydata/{}.dnds.info.csv".format(gene))
gwf.target_from_template("dNdSsummarybetween",
                         summary_files(dNdSfiles, "summaryfiles/divergent_genes.dnds.csv"))
gwf.target_from_template("pdistsummarybetween",
                         summary_files(pdistfiles, "summaryfiles/divergent_genes.pdist.csv"))
gwf.target_from_template("sitesummarybetween",
                         summary_files(sitefiles, "summaryfiles/divergent_genes.sites.csv"))
gwf.target_from_template("dNdSinfosummarybetween",
                         summary_files(dNdSinfofiles, "summaryfiles/divergent_genes.dnds.info.csv"))
pdistfiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    pdistfiles.append("summarydata_protein/{}.pdist.csv".format(gene))
sitefiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    sitefiles.append("summarydata_protein/{}.sites.csv".format(gene))
gwf.target_from_template("PEPpdistsummarybetween",
                         summary_files(pdistfiles, "summaryfiles/divergent_proteins.pdist.csv"))
gwf.target_from_template("PEPsitesummarybetween",
                         summary_files(sitefiles, "summaryfiles/divergent_proteins.sites.csv"))
############################################
## Sample random genes and run full pipeline
############################################
## LOAD random genes file: random_genes.txt
## Was created by:
## grep '>' occidentale/clover.cds.fa | python createRandomList.py random_genes.txt 30
# randomgenelist = readgenelist("random_genes.txt")
def SplitSampleFasta(filename, outdir, parts):
    inputs = [filename]
    outputs = []
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }
    fdir = "/".join(filename.split("/")[0:-1])
    for i in range(1, parts+1):
        if fdir=="":
            outputs.append(outdir+"/"+filename+".p"+str(i))
        else:
            outputs.append(fdir+"/"+outdir+"/"+filename.split("/")[-1]+".p"+str(i))
    spec = '''
    python splitfasta.py {} {} {}
    '''.format(filename, outdir, parts)
    return inputs, outputs, options, spec
parts = 100
gwf.target_from_template("SplitSamplegenes",
                         SplitSampleFasta("occidentale/clover.protein.fa", "parts", parts))
for i in range(1, parts+1):
    gwf.target_from_template("RepensTogenBlastPart"+str(i),
                         genblast("occidentale/parts/clover.protein.fa.p"+str(i), "repens/TrR.v5.To.fasta",
                                  "genblastresults/parts/Togenes.fasta.p"+str(i), "../occidentale/parts/TrTo.p"+str(i)))
    gwf.target_from_template("PalgenBlastPart"+str(i),
                             genblast("occidentale/parts/clover.protein.fa.p"+str(i), "pallescens/final_pallescens.fa",
                                      "genblastresults/parts/Palgenes.fasta.p"+str(i), "../occidentale/parts/Tp.p"+str(i)))
    gwf.target_from_template("PaltranslatePart"+str(i),
                             translateSequences("genblastresults/parts/Palgenes.fasta.p"+str(i), "genblastresults/parts/Palgenes.prot.fasta.p"+str(i)))
    gwf.target_from_template("RepensTpgenBlastPart"+str(i),
                             genblast("genblastresults/parts/Palgenes.prot.fasta.p"+str(i), "repens/TrR.v5.Tp.fasta",
                                      "genblastresults/parts/Tpgenes.fasta.p"+str(i), "../genblastresults/parts/TrTp.p"))
def join_files(sequences, outfile):
    inputs = sequences
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '00:20:00'
    }
    spec = "cat"
    for sequence in sequences:
        spec += " {} ".format(sequence)
    spec += "> {}".format(outfile)
    return inputs, outputs, options, spec
TrTogenes = []
for i in range(1, parts+1):
    TrTogenes.append("genblastresults/parts/Togenes.fasta.p"+str(i))
Palgenes = []
for i in range(1, parts+1):
    Palgenes.append("genblastresults/parts/Palgenes.fasta.p"+str(i))
TrTpgenes = []
for i in range(1, parts+1):
    TrTpgenes.append("genblastresults/parts/Tpgenes.fasta.p"+str(i))
gwf.target_from_template("AllRepensTo",
                         join_files(TrTogenes, "genblastresults/Togenes.all.fasta"))
gwf.target_from_template("AllRepensPal",
                         join_files(Palgenes, "genblastresults/Palgenes.all.fasta"))
gwf.target_from_template("AllRepensTp",
                         join_files(TrTpgenes, "genblastresults/Tpgenes.all.fasta"))
def efficient_workflow(Togenes, TrTogenes, Tpgenes, Palgenes, gene):
    inputs = [Togenes, TrTogenes, Tpgenes, Palgenes]
    outputs = []
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }
    spec = '''
    source activate python2
    echo "Filtering the gene from all of the respective species"
    python FASTafilter.py '_{gene}' {Togenes} > /scratch/$SLURM_JOBID/To00.gene.fasta
    python FASTafilter.py '_{gene}' {TrTogenes} > /scratch/$SLURM_JOBID/TrTo.gene.fasta
    python FASTafilter.py '_{gene}' {Tpgenes} > /scratch/$SLURM_JOBID/TrTp.gene.fasta
    python FASTafilter.py '_{gene}' {Palgenes} > /scratch/$SLURM_JOBID/Tp00.gene.fasta
    echo "Joining the gene files"
    python joingenes.py /scratch/$SLURM_JOBID/To00.gene.fasta /scratch/$SLURM_JOBID/TrTo.gene.fasta /scratch/$SLURM_JOBID/TrTp.gene.fasta /scratch/$SLURM_JOBID/Tp00.gene.fasta > /scratch/$SLURM_JOBID/joinedgenes.fasta
    cat /scratch/$SLURM_JOBID/joinedgenes.fasta
    ## Join the genes and translate
    python joingenesandtranslate.py  /scratch/$SLURM_JOBID/To00.gene.fasta /scratch/$SLURM_JOBID/TrTo.gene.fasta /scratch/$SLURM_JOBID/TrTp.gene.fasta /scratch/$SLURM_JOBID/Tp00.gene.fasta > /scratch/$SLURM_JOBID/joinedgenes.prot.fasta
    cat /scratch/$SLURM_JOBID/joinedgenes.prot.fasta
    ## Change environment
    source activate prank
    ## Align genes
    echo "Aligning the genes"
    prank -d=/scratch/$SLURM_JOBID/joinedgenes.fasta -o=/scratch/$SLURM_JOBID/alignedgenes -codon +F
    prank -d=/scratch/$SLURM_JOBID/joinedgenes.prot.fasta -o=/scratch/$SLURM_JOBID/alignedgenes.prot -protein
    cat /scratch/$SLURM_JOBID/alignedgenes.best.fas
    cat /scratch/$SLURM_JOBID/alignedgenes.prot.best.fas
    source activate python2
    ## Calculate measure on the alignments
    python dNdS.py /scratch/$SLURM_JOBID/alignedgenes.best.fas /scratch/$SLURM_JOBID/dnds.csv /scratch/$SLURM_JOBID/pdist.csv /scratch/$SLURM_JOBID/sites.csv /scratch/$SLURM_JOBID/dnds.sites.csv
    python comparePeptides.py /scratch/$SLURM_JOBID/alignedgenes.prot.best.fas /scratch/$SLURM_JOBID/pdist.prot.csv /scratch/$SLURM_JOBID/sites.prot.csv
    cat /scratch/$SLURM_JOBID/dnds.csv
    cat /scratch/$SLURM_JOBID/pdist.csv
    cat /scratch/$SLURM_JOBID/sites.csv
    cat /scratch/$SLURM_JOBID/dnds.sites.csv
    cat /scratch/$SLURM_JOBID/pdist.prot.csv
    cat /scratch/$SLURM_JOBID/sites.prot.csv
    python printcsv.py /scratch/$SLURM_JOBID/dnds.csv "{gene}" all_summaryfiles/all.genes.dnds.csv
    python printcsv.py /scratch/$SLURM_JOBID/pdist.csv "{gene}" all_summaryfiles/all.genes.pdist.csv
    python printcsv.py /scratch/$SLURM_JOBID/sites.csv "{gene}" all_summaryfiles/all.genes.sites.csv
    python printcsv.py /scratch/$SLURM_JOBID/dnds.sites.csv "{gene}" all_summaryfiles/all.genes.dnds.sites.csv
    python printcsv.py /scratch/$SLURM_JOBID/pdist.prot.csv "{gene}" all_summaryfiles/all.genes.pdist.prot.csv
    python printcsv.py /scratch/$SLURM_JOBID/sites.prot.csv "{gene}" all_summaryfiles/all.genes.sites.prot.csv
    '''.format(Togenes=Togenes, TrTogenes=TrTogenes, Tpgenes=Tpgenes, Palgenes=Palgenes, gene=gene)
    return inputs, outputs, options, spec
fullgenelist = readgenelist("allgenes.txt")
#fullgenelist = occi_gl
#fullgenelist = ["56g36650.1", "209g75490.3", "172g65010.1", "912g20040.2", "726g61350.1", "558g23290.1"]
for gene in fullgenelist:
    gwf.target_from_template("Gene"+gene,
                             efficient_workflow("occidentale/clover.cds.fa",
                                                "genblastresults/Togenes.all.fasta",
                                                "genblastresults/Tpgenes.all.fasta",
                                                "genblastresults/Palgenes.all.fasta", gene))
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

Reads a fasta file and cleans the contents of any characters that should not be there.

``` python
import sys
def read_and_fix(filename):
    f = open(filename)
    inseq = False
    alphabet = set([chr(i) for i in xrange(65, 91)])
    while True:
        c = f.read(1)
        if not c: break
        if inseq:
            if c==">":
                inseq = False
                sys.stdout.write("\n")
                sys.stdout.write(c)
            if c.upper() not in alphabet: continue
        if c=="\n":
            inseq = True
        sys.stdout.write(c)
    print
if __name__=="__main__":
    read_and_fix(sys.argv[1])
```

FASTafilter
-----------

Fasta fitler script, which finds a given sequence using the sequence name. First using grep, then jump to the file location to readout the sequence.

``` python
import subprocess
from sys import argv, exit
def findlocation(filename, search_term):
    out = subprocess.check_output("grep -b '{}' {}".format(search_term, filename), shell=True)
    matches = out.split("\n")[:-1]
    if len(matches)>1:
        return [int(match.split(":")[0]) for match in matches]
    return int(out.split(":")[0])
def collectsequence(filename, location):
    f = open(filename)
    f.seek(location)
    filebegun = False
    s = ""
    for line in f:
        if filebegun==False:
            if ">" in line:
                filebegun = True
                s += line
        else:
            if ">" in line: break
            s += line
    return s
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
            if ">" in line: break
            print line,
if __name__=="__main__":
    search_term, filename = argv[1:3]
    location = findlocation(filename, search_term)
    if type(location)==list:
        seqs = []
        for loc in location:
            #seqs.append([collectsequence(filename, loc), loc])
            filtersequence(filename, loc)
        #location = sorted(seqs, key=lambda x: len(x[0]))[-1][1]
    else:
        filtersequence(filename, location)
```

Joins a list of given genes, and cleans the sequence names.

``` python
import sys
from CorrectORF import readfasta
import os
from stat import S_ISFIFO
def collect_genes(genes):
    gene_files = []
    for gene in genes:
        sequence = readfasta(gene).values()
        if sequence:
            gene_files.append([gene, sequence[0].upper()])
    return gene_files
#def correct_genes(gene_files, reference):
#    for i in range(1, len(gene_files)):
#        gene_files[i][1] = clean_seq(gene_files[i][1].upper(), reference)
#    return gene_files
def write_out_fasta(genes):
    for gene in genes:
        print ">"+gene[0].split(".")[0].split("/")[-1]
        print gene[1]
if __name__=="__main__":
    if S_ISFIFO(os.fstat(0).st_mode):
        reference = "".join(sys.stdin.readlines()[1:]).replace("\n","")
        genes = sys.argv[1:]
        gene_files = collect_genes(genes)
        #gene_files = correct_genes(gene_files, reference)
        write_out_fasta(gene_files)
    else:
        genes = sys.argv[1:]
        gene_files = collect_genes(genes)
        #gene_files = correct_genes(gene_files)
        write_out_fasta(gene_files)
## COLLECT THE SEQUENCES
## REFERENCE IS ALWAYS THE FIRST SEQUENCE
## CORRECT ALL OF THE SEQUENCES ACCORDING TO THE REFERENCE
## WRITE OUT THE GENES.
```

Joins a list of given genes, and cleans the sequence names. It also translate the sequence into protein.

``` python
import sys
from CorrectORF import *
import os
from stat import S_ISFIFO
def collect_genes(genes):
    gene_files = []
    for gene in genes:
        sequence = readfasta(gene).values()
        if sequence:
            gene_files.append([gene, translate(sequence[0].upper())])
    return gene_files
def correct_genes(gene_files, reference):
    for i in range(1, len(gene_files)):
        gene_files[i][1] = clean_seq(gene_files[i][1].upper(), reference)
    return gene_files
def write_out_fasta(genes):
    for gene in genes:
        print ">"+gene[0].split(".")[0].split("/")[-1]
        print gene[1]
if __name__=="__main__":
    if S_ISFIFO(os.fstat(0).st_mode):
        reference = "".join(sys.stdin.readlines()[1:]).replace("\n","")
        genes = sys.argv[1:]
        gene_files = collect_genes(genes)
        #gene_files = correct_genes(gene_files, reference)
        write_out_fasta(gene_files)
    else:
        genes = sys.argv[1:]
        gene_files = collect_genes(genes)
        #gene_files = correct_genes(gene_files)
        write_out_fasta(gene_files)
## COLLECT THE SEQUENCES
## REFERENCE IS ALWAYS THE FIRST SEQUENCE
## CORRECT ALL OF THE SEQUENCES ACCORDING TO THE REFERENCE
## WRITE OUT THE GENES.
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
    if pathways==None:
        return 0, 1
    if type(pathways[0])==int:
        pathways = [pathways]
    c = 0
    for pathway in pathways:
        if pathway == None: continue
        s += pathway[0]
        n += pathway[1]
        c += 1
    if c==0:
        return 0, 1
    return s/float(c), n/float(c)
def jukes_cantor(D):
    if D>=0.745: return 27
    return -3/float(4)*log(1-4/float(3)*D)
def countN_S(seq):
    seq_n = len(seq)
    print "Sequence length:", seq_n
    N, S = 0, 0
    otherbases = {'A': ['T', 'G', 'C'],
                  'T': ['A', 'C', 'G'],
                  'G': ['A', 'T', 'C'],
                  'C': ['A', 'T', 'G']}
    for i in range(0, seq_n, 3):
        codon = seq[i:i+3]
        n, s = 0, 0
        for j in range(0,3):
            aa = codon_map.get(codon, "-")
            if aa=="-": continue
            tn, ts = 0, 0
            for k in otherbases.get(codon[j], []):
                ncodon = codon[:j]+k+codon[j+1:]
                if aa==codon_map[ncodon]:
                    ts += 1
                else:
                    tn += 1
            n += tn/float(3)
            s += ts/float(3)
        N += n
        S += s
    return N, S
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
            #print codon1, codon2
            sdnd = count_pathways(codon1, codon2)
            #print sdnd
            Sd += sdnd[0]
            Nd += sdnd[1]
    N1, S1 = countN_S(seq1)
    N2, S2 = countN_S(seq2)
    N, S = (N1+N2)/float(2), (S1+S2)/float(2)
    if S==0 and N==0:
        return "NaN", (N, S)
    print "S & N", S, N
    print "Sd & Nd", Sd, Nd
    pS = Sd/float(S)
    pN = Nd/float(N)
    print "pS & pN", pS, pN
    #print "pN/pS ratio", pN/pS
    dS = jukes_cantor(pS)
    dN = jukes_cantor(pN)
    print "dS & dN", dS, dN
    if dS==0:
        return "NaN", (N, S)
    print "dN/dS ratio", dN/dS
    return dN/dS, (Nd, Sd, N, S)
def heterozygote_sites(seq):
    s = seq.replace("-", "")
    for char in ["A", "T", "G", "C", "N"]:
        s = s.replace(char, "")
    print s
    return len(s)
def pdistance(seq1, seq2):
    n, d = 0, 0
    for base1, base2 in zip(seq1.upper(), seq2.upper()):
        if base1=="-" or base2=="-": continue
        if base1!=base2:
            d += 1
        n += 1
    if n==0: return "NaN", (0,0)
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
    return keys, new_keys
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
                print(matrix[key][key2])
                f.write(", ({}".format(matrix[key][key2][0]))
                for item in matrix[key][key2][1:]:
                    f.write(";{}".format(item))
                f.write(")")
            else:
                f.write(", {}".format(matrix[key][key2]))
        f.write("\n")
    f.close()
if __name__=="__main__":
    sequences = readfasta(argv[1])
    done = []
    print sequences
    keys, samples = cleanKeyNames(sequences.keys())
    for key, sample in zip(keys, samples):
        sequences[sample] = sequences[key]
        #del sequences[key]
    print sequences
    for key in keys:
        if key not in samples:
            del sequences[key]
    print sequences
    dNdS = {}
    pdist = {}
    sites = {}
    dNdSinfo = {}
    for seq1 in sequences:
        if seq1 not in dNdS:
            dNdS[seq1] = {}
            pdist[seq1] = {}
            sites[seq1] = {}
            dNdSinfo[seq1] = {}
        for seq2 in sequences:
            dNdS[seq1][seq2], dNdSinfo[seq1][seq2] = calculate_dn_ds(sequences[seq1].upper(), sequences[seq2].upper())
            pdist[seq1][seq2], sites[seq1][seq2] = pdistance(sequences[seq1].upper(), sequences[seq2].upper())
    #prettyprint(dNdS)
    #prettyprint(pdist)
    write_out_matrix(argv[2], dNdS)
    write_out_matrix(argv[3], pdist)
    write_out_matrix(argv[4], sites)
    write_out_matrix(argv[5], dNdSinfo)
```

Joins a list of given genes, and cleans the sequence names.

``` python
from dNdS import *
if __name__=="__main__":
    sequences = readfasta(argv[1])
    done = []
    keys, samples = cleanKeyNames(sequences.keys())
    for key, sample in zip(sequences.keys(), samples):
        sequences[sample] = sequences[key]
        if key not in samples:
            del sequences[key]
    pdist = {}
    sites = {}
    for seq1 in sequences:
        if seq1 not in pdist:
            pdist[seq1] = {}
            sites[seq1] = {}
        for seq2 in sequences:
            pdist[seq1][seq2], sites[seq1][seq2] = pdistance(sequences[seq1], sequences[seq2])
    #prettyprint(dNdS)
    #prettyprint(pdist)
    write_out_matrix(argv[2], pdist)
    write_out_matrix(argv[3], sites)
```

Joins summary files returned from dNdS into one summary file.

``` python
from sys import argv
def read_csv_files(filenames):
    full_data = {}
    full_data["ID"] = []
    for filename in filenames:
        f = open(filename)
        for item in full_data:
            full_data[item].append("-")
        full_data["ID"][-1] = filename.split("/")[-1]
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
                full_data[item] = ["-" for i in range(len(full_data["ID"]))]
            full_data[item][-1] = row[item]
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

Gene annotation pipeline
========================

Gene annotation pipeline which uses BRAKER2 to perform the annotation. Comparing the BRAKER2 annotation based on AUGUSTUS using RNAseq evidence and Medicago geneset, to previous annotations. The previous annotations where a Medicago based annotation that used GMAP to map medicago genes, and the first annotation. The main workflow built using gwf. Uses local installations of BRAKER2, BUSCO and gmap. The final output is a in the form of a .gtf file which has been curated based on evidence from the other gene annotations we have, and a BUSCO report.

``` python
from gwf import Workflow
gwf = Workflow()
def extract_sequence(reference, sequence_name, output):
    """
    Extract sequence takes in a <reference>, <sequence_name> and <output>
    The <reference> is a fasta file, which must be indexed. (have a .fai) file
    It looks up in the .fai file to find the <sequence_name>.
    If the name is present it reads the sequences from the <reference>,
    and saves it as <output>
    """
    inputs = [reference]
    outputs = [output]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }
    directory = "/".join(output.split("/")[:-1])
    spec = '''
    mkdir -p {dirc}
    python extract_sequence.py {ref} {seq} > {out}
    '''.format(ref=reference, seq=sequence_name, dirc=directory, out=output)
    return inputs, outputs, options, spec
def repeatmasking(reference, output_reference, directory):
    """
    """
    inputs = [reference]
    outputs = [output_reference]
    options = {
        'cores': 8,
        'memory': '24g',
        'account': 'NChain',
        'walltime': '08:00:00'
    }
    spec = '''
    ../../ANNOTATION/RepeatMasker/bin/RepeatMasker -xsmall -pa {cores} -dir {dirc} -species arabidopsis {ref}
    '''.format(ref=reference, cores=options['cores'], dirc=directory)
    return inputs, outputs, options, spec
def fasta_index(reference):
    inputs = [reference]
    outputs = [reference+".fai"]
    options = {
        'cores': 1,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }
    spec = '''
    source /com/extra/samtools/1.6.0/load.sh
    samtools faidx {ref}
    '''.format(ref=reference)
    return inputs, outputs, options, spec
gwf.target_from_template("RepensMasking",
                         repeatmasking("repens/TrR.v5.fasta",
                                       "repens/TrR.v5.fasta.masked",
                                       "repens"))
gwf.target_from_template("OccidentaleMasking",
                         repeatmasking("occidentale/To.fasta",
                                       "occidentale/To.fasta.masked",
                                       "occidentale"))
gwf.target_from_template("PallescensMasking",
                         repeatmasking("pallescens/Tp.fasta",
                                       "pallescens/Tp.fasta.masked",
                                       "pallescens"))
gwf.target_from_template("RepensMaskindex",
                         fasta_index("repens/TrR.v5.fasta.masked"))
gwf.target_from_template("OccidentaleMaskindex",
                         fasta_index("occidentale/To.fasta.masked"))
gwf.target_from_template("PallescensMaskindex",
                         fasta_index("pallescens/Tp.fasta.masked"))
def get_sequence_names(reference):
    f = open(reference+".fai")
    sequence_names = []
    for line in f:
        sequence_names.append(line.split("\t")[0])
    return sequence_names
white_clover_chromosomes = get_sequence_names("repens/TrR.v5.fasta")
occidentale_chromosomes = get_sequence_names("occidentale/To.fasta")
pallescens_chromosomes = get_sequence_names("pallescens/Tp.fasta")
n = 10
print("White clover, chromosomes:", len(white_clover_chromosomes))
print("Occidentale, chromosomes:", len(occidentale_chromosomes))
print("Pallescens, chromosomes:", len(pallescens_chromosomes))
for chromosome in white_clover_chromosomes:
    gwf.target_from_template("Repens"+chromosome+"Extract",
                             extract_sequence("repens/TrR.v5.fasta.masked",
                                              chromosome,
                                              "runs/repens/"+chromosome+"/"+chromosome+".fa"))
#for chromosome in occidentale_chromosomes[:n]:
#     gwf.target_from_template("Occidentale"+chromosome+"Extract",
#                              extract_sequence("occidentale/To.fasta.masked",
#                                               chromosome,
#                                               "runs/occidentale/"+chromosome+"/"+chromosome+".fa"))
#for chromosome in pallescens_chromosomes[:n]:
#     gwf.target_from_template("Pallescens"+chromosome+"Extract",
#                              extract_sequence("pallescens/final_pallescens.fa.masked",
#                                               chromosome,
#                                               "runs/pallescens/"+chromosome+"/"+chromosome+".fa"))
gwf.target("OccidentaleSequence", inputs=["occidentale/To.fasta.masked"],
           outputs=["runs/occidentale/occidentalefull/To.fasta"]) << """
mkdir -p runs/occidentale/occidentalefull
cp occidentale/To.fasta.masked runs/occidentale/occidentalefull/To.fasta
"""
gwf.target("PallesecensSequence", inputs=["pallescens/Tp.fasta.masked"],
           outputs=["runs/pallescens/pallescensfull/Tp.fasta"]) << """
mkdir -p runs/pallescens/pallescensfull
cp pallescens/Tp.fasta.masked runs/pallescens/pallescensfull/Tp.fasta
"""
def filter_bam_file(bamfile, chromosome, outfile):
    """
    filter_bam_file uses samtools to read a <bamfile> and read only
    the reads that are mapped to <chromosome>.
    It saves the filtered reads into <outfile>.
    """
    inputs = [bamfile]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }
    directory = "/".join(outfile.split("/")[:-1])
    spec = '''
    source /com/extra/samtools/1.6.0/load.sh
    mkdir -p {dirc}
    samtools view -b {infile} {chrom} > {out}
    '''.format(infile=bamfile, chrom=chromosome, out=outfile, dirc=directory)
    return inputs, outputs, options, spec
for chromosome in white_clover_chromosomes:
    gwf.target_from_template("Repens"+chromosome+"BamExtract",
                             filter_bam_file("RNAdata/repens/pooled.reads.bam",
                                             chromosome,
                                             "runs/repens/"+chromosome+"/"+chromosome+".bam"))
# for chromosome in occidentale_chromosomes[:n]:
#      gwf.target_from_template("Occidentale"+chromosome+"BamExtract",
#                               filter_bam_file("RNAdata/occidentale/pooled_reads.bam",
#                                               chromosome,
#                                               "runs/occidentale/"+chromosome+"/"+chromosome+".bam"))
# for chromosome in pallescens_chromosomes[:n]:
#      gwf.target_from_template("Pallescens"+chromosome+"BamExtract",
#                               filter_bam_file("RNAdata/pallescens/pooled_reads.bam",
#                                               chromosome,
#                                               "runs/pallescens/"+chromosome+"/"+chromosome+".bam"))
gwf.target("OccidentaleReads", inputs=["RNAdata/occidentale/pooled_reads.bam"],
           outputs=["runs/occidentale/occidentalefull/pooled_reads.bam"]) << """
mkdir -p runs/occidentale/occidentalefull
cp RNAdata/occidentale/pooled_reads.bam runs/occidentale/occidentalefull/pooled_reads.bam
"""
gwf.target("PallesecensReads", inputs=["RNAdata/pallescens/pooled_reads.bam"],
           outputs=["runs/pallescens/pallescensfull/pooled_reads.bam"]) << """
mkdir -p runs/pallescens/pallescensfull
cp RNAdata/pallescens/pooled_reads.bam runs/pallescens/pallescensfull/pooled_reads.bam
"""
def BRAKER_gene_annotation(reference, bamfile, proteinfile, directory, species, outfile, time='12:00:00', cores=4):
    """
    Gene annoation using BRAKER.
    """
    inputs = [reference, bamfile, proteinfile]
    outputs = [outfile]
    options = {
        'cores': cores,
        'memory': '16g',
        'account': 'NChain',
        'walltime': time
    }
    levels = len(directory.split("/"))+1
    bamfile = bamfile.split("/")[-1]
    reference = reference.split("/")[-1]
    proteinfile = "../"*(levels-1)+proteinfile
    
    spec = '''
    cd {directory}
    
    source activate BRAKER
    export GENEMARK_PATH=/faststorage/project/NChain/WHITE_CLOVER/BRAKER_ANNOTATION/gm_et_linux_64/gmes_petap
    export ALIGNMENT_TOOL_PATH=/faststorage/project/NChain/WHITE_CLOVER/BRAKER_ANNOTATION/gth-1.7.1-Linux_x86_64-64bit/bin
    export AUGUSTUS_BIN_PATH=/home/marnit/anaconda3/envs/BRAKER/bin/
    export AUGUSTUS_CONFIG_PATH=/home/marnit/anaconda3/envs/BRAKER/config/
    export AUGUSTUS_SCRIPTS_PATH=/home/marnit/anaconda3/envs/BRAKER/scripts/
    export PERL5LIB=/home/marnit/anaconda3/envs/BRAKER/lib/site_perl/5.26.2
    export PERL5LIB="/home/marnit/anaconda3/envs/BRAKER/lib/site_perl/5.26.2/x86_64-linux-thread-multi:$PERL5LIB"
    export PATH="/home/marnit/anaconda3/envs/BRAKER/bin/:$PATH"
    {calldir}BRAKER_v2.1.0/braker.pl --genome={ref} --prot_seq={prot} --prg=gth --bam={bam} --cores={cores} --species={species} --softmasking --useexisting
    '''.format(directory=directory, calldir="../"*levels, species=species,
               cores=options['cores'], ref=reference, bam=bamfile, prot=proteinfile)
    return inputs, outputs, options, spec
for chromosome in white_clover_chromosomes:
    gwf.target_from_template("Repens"+chromosome+"Annotation",
                             BRAKER_gene_annotation("runs/repens/"+chromosome+"/"+chromosome+".fa",
                                                    "runs/repens/"+chromosome+"/"+chromosome+".bam",
                                                    "MedicagoProtein/Mt.proteins.fa",
                                                    "runs/repens/"+chromosome,
                                                    "repens"+chromosome,
                                                    "runs/repens/"+chromosome+"/braker/repens"+chromosome+"/augustus.hints.gtf"))
# for chromosome in occidentale_chromosomes[:n]:
#     gwf.target_from_template("Occidentale"+chromosome+"Annotation",
#                              BRAKER_gene_annotation("runs/occidentale/"+chromosome+"/"+chromosome+".fa",
#                                                     "runs/occidentale/"+chromosome+"/"+chromosome+".bam",
#                                                     "MedicagoProtein/Mt.proteins.fa",
#                                                     "runs/occidentale/"+chromosome,
#                                                     "occidentale"+chromosome,
#                                                     "runs/occidentale/"+chromosome+"/braker/occidentale"+chromosome+"/augustus.hints.gtf"))
# for chromosome in pallescens_chromosomes[:n]:
#     gwf.target_from_template("Pallescens"+chromosome+"Annotation",
#                              BRAKER_gene_annotation("runs/pallescens/"+chromosome+"/"+chromosome+".fa",
#                                                     "runs/pallescens/"+chromosome+"/"+chromosome+".bam",
#                                                     "MedicagoProtein/Mt.proteins.fa",
#                                                     "runs/pallescens/"+chromosome,
#                                                     "pallescens"+chromosome,
#                                                     "runs/pallescens/"+chromosome+"/braker/pallescens"+chromosome+"/augustus.hints.gtf"))
gwf.target_from_template("OccidentaleAnnotation",
                         BRAKER_gene_annotation("runs/occidentale/occidentalefull/To.fasta",
                                                "runs/occidentale/occidentalefull/pooled_reads.bam",
                                                "MedicagoProtein/Mt.proteins.fa",
                                                "runs/occidentale/occidentalefull",
                                                "occidentale",
                                                "runs/occidentale/occidentalefull/braker/occidentale/augustus.hints.gtf",
                                                "48:00:00",
                                                8))
gwf.target_from_template("PallescensAnnotation",
                         BRAKER_gene_annotation("runs/pallescens/pallescensfull/Tp.fasta",
                                                "runs/pallescens/pallescensfull/pooled_reads.bam",
                                                "MedicagoProtein/Mt.proteins.fa",
                                                "runs/pallescens/pallescensfull",
                                                "pallescens",
                                                "runs/pallescens/pallescensfull/braker/pallescens/augustus.hints.gtf",
                                                "48:00:00",
                                                8))
def join_hints(hints, output):
    """ """
    inputs = hints
    outputs = [output]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }
    spec = "cat"
    for hint in hints:
        spec += " {}".format(hint)
    spec += " > {}".format(output)
    return inputs, outputs, options, spec
repens_hints = []
for chromosome in white_clover_chromosomes[1:]:
    repens_hints.append("runs/repens/"+chromosome+"/braker/repens"+chromosome+"/augustus.hints.gtf")
gwf.target_from_template("JoinRepensAnnotation",
                         join_hints(repens_hints,
                                    "runs/repens/repens.annotation.gtf"))
# occidentale_hints = []
# for chromosome in occidentale_chromosomes[:n]:
#      occidentale_hints.append("runs/occidentale/"+chromosome+"/braker/occidentale"+chromosome+"/augustus.hints.gtf")
# gwf.target_from_template("JoinOccidentaleAnnoation",
#                          join_hints(occidentale_hints,
#                                     "runs/occidentale/occidentale.annotation.gtf"))
# pallescens_hints = []
# for chromosome in pallescens_chromosomes[:n]:
#      pallescens_hints.append("runs/pallescens/"+chromosome+"/braker/pallescens"+chromosome+"/augustus.hints.gtf")
# gwf.target_from_template("JoinPallescensAnnoation",
#                           join_hints(pallescens_hints,
#                                      "runs/pallescens/pallescens.annotation.gtf"))
gwf.target("OccidentaleFinal", inputs=["runs/occidentale/occidentalefull/braker/occidentale/augustus.hints.gtf"],
           outputs=["runs/occidentale/occidentale.annotation.gtf"]) << """
cp runs/occidentale/occidentalefull/braker/occidentale/augustus.hints.gtf runs/occidentale/occidentale.annotation.gtf
"""
gwf.target("PallesecensFinal", inputs=["runs/pallescens/pallescensfull/braker/pallescens/augustus.hints.gtf"],
           outputs=["runs/pallescens/pallescens.annotation.gtf"]) << """
cp runs/pallescens/pallescensfull/braker/pallescens/augustus.hints.gtf runs/pallescens/pallescens.annotation.gtf
"""
def GMAP_prepare(reference, dbname):
    """ """
    inputs = [reference]
    outputs = [dbname]
    options = {
        'cores': 1,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }
    directory = dbname.split("/")[0]
    spec = '''
    source /com/extra/gmap/2017-08-15/load.sh
    mkdir -p {dirc}
    gmap_build -d {dbname} {ref} -D $PWD/{dirc}
    '''.format(dbname=dbname.split("/")[-1], ref=reference, dirc=directory)
    return inputs, outputs, options, spec
gwf.target_from_template("OccidentaleDB",
                         GMAP_prepare("occidentale/final.assembly.fa", "gmap/occidentale.gmap"))
gwf.target_from_template("PallescensDB",
                         GMAP_prepare("pallescens/final_pallescens.fa", "gmap/pallescens.gmap"))
def GMAP_mapping(dbname, genes, gfffile):
    """ """
    inputs = [dbname, genes]
    outputs = [gfffile]
    options = {
        'cores': 8,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '24:00:00'
    }
    directory = dbname.split("/")[0]
    spec = '''
    source /com/extra/gmap/2017-08-15/load.sh
    gmap -d {} -D $PWD/{} -t 8 -f 2 {} > {}
    '''.format(dbname.split("/")[1], directory, genes, gfffile, gfffile, gfffile)
    return inputs, outputs, options, spec
gwf.target_from_template("OccidentaleMedicago",
                         GMAP_mapping("gmap/occidentale.gmap", "MedicagoProtein/Mt4.0v2_Genes.fasta",
                                      "gmap/occidentale.Mt.gff"))
gwf.target_from_template("PallescensMedicago",
                         GMAP_mapping("gmap/pallescens.gmap", "MedicagoProtein/Mt4.0v2_Genes.fasta",
                                      "gmap/pallescens.Mt.gff"))
def gff3plsorting(gff_file, outgff_file):
    """ """
    inputs = [gff_file]
    outputs = [outgff_file]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }
    spec = '''
    /home/marnit/bin/gff3sort.pl --chr_order natural {infile} > {outfile}
    '''.format(infile=gff_file, outfile=outgff_file)
    return inputs, outputs, options, spec
gwf.target_from_template("OccidentaleMtSort",
                         gff3plsorting("gmap/occidentale.Mt.gff", "gff_files/occidentale.Mt.gff"))
gwf.target_from_template("PallescensMtSort",
                         gff3plsorting("gmap/pallescens.Mt.gff", "gff_files/pallescens.Mt.gff"))
def copy_files(infiles, outfiles):
    """ """
    inputs = infiles
    outputs = outfiles
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }
    spec = ''
    for inf, outf in zip(infiles, outfiles):
        spec += "cp {} {}\n".format(inf,outf)
    return inputs, outputs, options, spec
    
gwf.target_from_template("MovingFinalGFFfiles",
                         copy_files(["runs/repens/repens.annotation.gtf",
                                     "runs/occidentale/occidentale.annotation.gtf",
                                     "runs/pallescens/pallescens.annotation.gtf"],
                                    ["gff_files/repens.braker.gtf",
                                     "gff_files/occidentale.braker.gtf",
                                     "gff_files/pallescens.braker.gtf"]))
def ANNOTATIONcleaning(main_gff_file, additional_gff_files, output, confidence_level, binsize=25000):
    """ """
    inputs = [main_gff_file]+additional_gff_files
    outputs = [output]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }
    spec = '''
    source activate python2
python ANNOTATIONmerger.py -h
python ANNOTATIONcleaner.py -o {outfile} -c {conflvl} -b {binsize} {main_gff}'''.format(outfile=output, conflvl=confidence_level,
                                                                                        binsize=binsize, main_gff=main_gff_file)
    for agf in additional_gff_files:
        spec += " {}".format(agf)
    return inputs, outputs, options, spec
gwf.target_from_template("RepensCleaning",
                         ANNOTATIONcleaning("gff_files/repens.braker.gtf",
                                            ["gff_files/TrR.Mt.gff", "gff_files/TrR.v5.sorted.gff"],
                                            "gff_files/repens.final.gtf",
                                            0.25))
gwf.target_from_template("OccidentaleCleaning",
                         ANNOTATIONcleaning("gff_files/occidentale.braker.gtf",
                                            ["gff_files/To.v5.gff", "gff_files/occidentale.Mt.gff"],
                                            "gff_files/occidentale.final.gtf",
                                            0.25))
gwf.target_from_template("PallescensCleaning",
                         ANNOTATIONcleaning("gff_files/pallescens.braker.gtf",
                                            ["gff_files/Tp.v4.gff", "gff_files/pallescens.Mt.gff"],
                                            "gff_files/pallescens.final.gtf",
                                            0.25))
def BUSCO(ingff, reference, proteinfile, dbname, report):
    """ """
    inputs = [ingff, reference]
    outputs = [proteinfile, "busco/"+dbname]
    options = {
        'cores': 4,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }
    proteinfilename = proteinfile.split("/")[-1]
    
    spec = '''
    source /com/extra/cufflinks/2.2.1/load.sh
    source activate busco
    gffread -y {protfile} -g {ref} {ingff}
    
    cd busco
    run_BUSCO.py -i {protfilename} -o {name} -l ../../busco/sample_data/embryophyta_odb9/ -m proteins -c 4 > {report}
    '''.format(protfile=proteinfile, protfilename=proteinfilename, ref=reference, ingff=ingff, name=dbname, report=report)
    return inputs, outputs, options, spec
gwf.target_from_template("RepensBUSCO",
                         BUSCO("gff_files/repens.final.gtf",
                               "repens/TrR.v5.fasta",
                               "busco/repens.proteins.fa",
                               "RepensFinal",
                               "RepensFinalReport.out"))
gwf.target_from_template("OccidentaleBUSCO",
                         BUSCO("gff_files/occidentale.final.gtf",
                               "occidentale/final.assembly.fa",
                               "busco/occidentale.proteins.fa",
                               "OccidentaleFinal",
                               "OccidentaleFinalReport.out"))
gwf.target_from_template("PallescensBUSCO",
                         BUSCO("gff_files/pallescens.final.gtf",
                               "pallescens/final_pallescens.fa",
                               "busco/pallescens.proteins.fa",
                               "PallescensFinal",
                               "PallescensFinalReport.out"))
```

Efficient sequence reader from a fasta file. Used to seperate the chromosomes into seperate files to speed up the pipeline.

``` python
from sys import argv, exit, stdout
def read_index(filename):
    f = open(filename)
    indices = {}
    for line in f:
        chrom, _, byteindex, _, _ = line.split("\t")
        indices[chrom] = int(byteindex)
    return indices
def read_sequence(filename, index, target):
    f = open(filename)
    f.seek(index-len(target)-2)
    while True:
        c = f.read(1)
        if c==">":
            stdout.write(">")
            break
    stdout.write(f.readline())
    for line in f:
        if line[0]==">": break
        stdout.write(line)
if __name__=="__main__":
    if len(argv)==3:
        reference, target = argv[1:]
    else:
        print("USAGE: python extract_sequence.py <reference> <sequence_name>")
        exit()
    try:
        indices = read_index(reference+".fai")
    except:
        print("No index file included or file does not exist")
        exit()
    if target in indices:
        read_sequence(reference, indices[target], target)
    else:
        print("Sequence name not found!")
        exit()
```

Loads all of the exons from the additional gff files in to a hash(dictionary) for exact matches. For inexact it builds an exon map with hash lookup based on a binsize, and then compares the overlap of all the exons inside the bin. Computes a confidence based on exact and inexact matches, weighing more on exact matches (all exact also have inexact matches).

``` python
from ANNOTATIONmerger import *
from math import floor
def overlap(exon1, exon2):
    _, start1, stop1, dirc1 = exon1
    _, start2, stop2, dirc2 = exon2
    dirc = dirc1==dirc2
    left = stop2>start1 and start2<=start1
    right = start2<stop1 and stop2>=stop1
    in1 = start1<=start2 and stop1>=stop2
    in2 = start2<=start1 and stop2>=stop1
    if dirc and (left or right or in1 or in2):
        return 1
    else:
        return 0
def overlap_confidence(exon1, exon2):
    _, start1, stop1, dirc1 = exon1
    _, start2, stop2, dirc2 = exon2
    if dirc1!=dirc2: return 0
    area = min(stop1, stop2)-max(start1, start2)
    if area<0: area = 0
    return area
if __name__=="__main__":
    _current_version = "0.2.7"
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m <main gff file> <gff files to filter by>\033[39m'''
    parser = OptionParser(usage)
    #parser.add_option('-R', type="string", nargs=1, dest="RNAdepth",
    #                  help="RNA read depth file. <Created by samtools depth -aa reads.bam>")
    #parser.add_option('-Q', type="int", nargs=1, dest="QueueSize", default=15,
    #                  help="Maximum Queue size. Ability to look back to previous transcripts. Larger Q increases runtime and space consumption. Also allows for better comparisons of iso-forms. Default is set to 15.")
    parser.add_option('-c', type="string", nargs=1, dest="Cutoff", default="0.5",
                      help="Confidence level cutoff. Default is set to 0.5.")
    parser.add_option('-b', type="string", nargs=1, dest="BinSize", default="25000",
                      help="bin size. Default is set to 25 kb (30000).")
    parser.add_option('-o', type="string", nargs=1, dest="Output", default="output.gtf",
                      help="Output file name. Default: output.gtf")
    pretty_print("======= ANNOTATION CLEANER (v{}) =======".format(_current_version), "orange")
    
    options, args = parser.parse_args()
    
    _cut_off = float(options.Cutoff)
    _outfile_name = options.Output
    _bin_size = int(options.BinSize)
    _queue_size = 1
    pretty_print("Cutoff level set to {} and Bin size set to {}".format(_cut_off, _bin_size), "red")
    pretty_print("Output will be saved into: {}".format(_outfile_name), "red")
    pretty_print("OPENING GFF/GTF FILES", "lightblue")
    gene_annotations = []
    filenames = []
    for arg in args:
        if ".gff" in arg or ".gtf" in arg:
            filenames.append(arg)
    main_gff = GFF(filenames[0])
    pretty_print("\tMain gff file: "+filenames[0])
            
    for filename in filenames[1:]:
        pretty_print("\tAdditional gff files: "+filename)
        gene_annotations.append(GFF(filename))
    pretty_print("READING EXONS FROM THE ADDITIONAL FILES", "lightblue")
    exon_hash = {}
    exon_map = {}
    display_counter = 1
    display_by = 25000
    for ga in gene_annotations:
        while not ga._empty:
            if (display_counter % display_by)==0:
                print("\t{} transcripts read through".format(display_counter))
            gene = ga.get_gene()
            if 'exon' in gene.elements:
                for exon in gene.elements['exon']:
                    string_representation = "%r" % (exon)
                    exon_hash[string_representation] = 0
                    chrom, start, stop, dirc = exon
                    
                    if chrom not in exon_map:
                        exon_map[chrom] = {}
                    if dirc not in exon_map[chrom]:
                        exon_map[chrom][dirc] = {}
                    sbin = int(floor(start/_bin_size))*_bin_size
                    stbin = int(floor(stop/_bin_size))*_bin_size
                    if sbin==stbin:
                        if sbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][sbin] = []
                        exon_map[chrom][dirc][sbin].append(exon)
                    else:
                        if sbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][sbin] = []
                        exon_map[chrom][dirc][sbin].append(exon)
                        if stbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][stbin] = []
                        exon_map[chrom][dirc][stbin].append(exon)
                    
            if 'cds' in gene.elements:
                for exon in gene.elements['cds']:
                    string_representation = "%r" % (exon)
                    exon_hash[string_representation] = 0
                    chrom, start, stop, dirc = exon
                    
                    if chrom not in exon_map:
                        exon_map[chrom] = {}
                    if dirc not in exon_map[chrom]:
                        exon_map[chrom][dirc] = {}
                    sbin = int(floor(start/_bin_size))*_bin_size
                    stbin = int(floor(stop/_bin_size))*_bin_size
                    if sbin==stbin:
                        if sbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][sbin] = []
                        exon_map[chrom][dirc][sbin].append(exon)
                    else:
                        if sbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][sbin] = []
                        exon_map[chrom][dirc][sbin].append(exon)
                        if stbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][stbin] = []
                        exon_map[chrom][dirc][stbin].append(exon)
            display_counter += 1
    pretty_print("FILTERING MAIN GFF FILE", "lightblue")
    outgff = GFFwriter(_outfile_name, _current_version)
    total_gene_count = 0
    passed_genes = 0
    passed_genes_sub1 = 0
    passed_genes_sub2 = 0
    display_counter = 1
    display_by = 25000
    subgenome1 = set(['chr1', 'chr2', 'chr3', 'chr4',
                      'chr5', 'chr6', 'chr7', 'chr8'])
    subgenome2 = set(['chr9', 'chr10', 'chr11', 'chr12',
                      'chr13', 'chr14', 'chr15', 'chr16'])
    
    while not main_gff._empty:
        if (display_counter % display_by)==0:
                print("\t{} transcripts read through".format(display_counter))
        total_gene_count += 1
        gene = main_gff.get_gene()
        largest_type = max((len(gene.get('cds', [])), 'cds'),
                           (len(gene.get('exon', [])), 'exon'),
                           key=lambda x: x[0])[1]
        exon_count = 0
        exon_overlap_count = 0
        exon_overlap_score = 0
        for exon in gene[largest_type]:
            if "%r" % (exon) in exon_hash:
                exon_count += 1
            chrom, start, stop, dirc = exon
            if chrom in exon_map and dirc in exon_map.get(chrom, {}):
                sbin = int(floor(start/_bin_size))*_bin_size
                stbin = int(floor(stop/_bin_size))*_bin_size
                if sbin==stbin:
                    for test_exon in exon_map[chrom][dirc].get(sbin, []):
                        if overlap(exon, test_exon):
                            exon_overlap_count += 1
                            break
                        #exon_overlap_score += overlap_confidence(exon, test_exon)
                else:
                    for test_exon in exon_map[chrom][dirc].get(sbin, []):
                        if overlap(exon, test_exon):
                            exon_overlap_count += 1
                            break
                    for test_exon in exon_map[chrom][dirc].get(stbin, []):
                        if overlap(exon, test_exon):
                            exon_overlap_count += 1
                            break
        
        confidence = (exon_count+exon_overlap_count)/float(2*len(gene[largest_type]))
        
        if confidence>=_cut_off:
            passed_genes += 1
            if gene.chromosome in subgenome1:
                passed_genes_sub1 += 1
            if gene.chromosome in subgenome2:
                passed_genes_sub2 += 1
            outgff.write_gene(gene, confidence, filenames[0])
        display_counter += 1
    pretty_print("FILE STATISTICS", "lightblue")
    pretty_print("\tTotal gene number: {}".format(total_gene_count))
    pretty_print("\tPassed genes: {}".format(passed_genes))
    pretty_print("\tPassed genes in first subgenome: {}".format(passed_genes_sub1))
    pretty_print("\tPassed genes in second subgenome: {}".format(passed_genes_sub2))
    pretty_print("\tTotal: {}".format(passed_genes_sub1+passed_genes_sub2))
```

First attempt at doing a merging of the multiple gff files. Currently only used for the classes and functions, inside on ANNOTATIONcleaner.py

``` python
#import sys
from optparse import OptionParser
#from matplotlib import pyplot
class Queue:
    def __init__(self, queue_size):
        self._list = list()
        self._max_size = queue_size
        self.full = False
    def add(self, value):
        if self.full:
            self._list.pop(0)
            self._list.append(value)
        else:
            self._list.append(value)
            if len(self._list)==self._max_size:
                self.full = True
    def get(self, element):
        if element>self._max_size: return None
        if len(self._list)<element: return None
        return self._list[-element]
class Gene:
    def __init__(self, information):
        self.elements = {}
        self.new_gene(information)
    def new_gene(self, information):
        chrom, _, element_type, start, stop, _, direction = information[:7]
        self.elements['transcript'] = [chrom, int(start), int(stop), direction]
        self.direction = direction
        self.start = int(start)
        self.stop = int(stop)
        self.chromosome = chrom
        if "_" in chrom:
            self.chrom_number = int(chrom.split("_")[-1])
        else:
            self.chrom_number = int(chrom.split("chr")[-1])
        self._lines = [information]
    def add_element(self, information):
        chrom, _, element_type, start, stop, _, direction = information[:7]
        element_type = element_type.lower()
        if self.direction==".":
            self.direction = direction
        #print self.__contains__([chrom, start, stop, direction])
        if [chrom, int(start), int(stop), direction] in self:
            self._lines.append(information)
            if element_type.lower() not in self.elements:
                self.elements[element_type.lower()] = list()
            self.elements[element_type.lower()].append([chrom, int(start), int(stop), direction])
    def __contains__(self, gene_info):
        if type(gene_info)==Gene:
            return gene_info.direction==self.direction and gene_info.chromosome==self.chromosome and (self.start>=gene_info.start and self.stop<=gene_info.stop)
        if type(gene_info)==list or type(gene_info)==tuple:
            if len(gene_info)==4:
                c, start, stop, d = gene_info
                return d==self.direction and c==self.chromosome and self.start<=start and self.stop>=stop
            else:
                print "WARNING: If list input must be [chromosome, start, stop, direction]"
                return None
    def overlaps(self, gene):
        chrom = gene.chromosome==self.chromosome
        dirc = gene.direction==self.direction
        left = gene.stop>self.start and gene.start<=self.start
        right = gene.start<self.stop and gene.stop>=self.stop
        return (chrom and dirc) and (left or right or (gene in self) or (self in gene))
    def __lt__(self, gene):
        if self.chromosome!=gene.chromosome:
            return self.chrom_number<gene.chrom_number
        return self.stop<gene.stop
    def __le__(self, gene):
        if self.chromosome!=gene.chromosome:
            return self.chrom_number<=gene.chrom_number
        return self.stop<=gene.stop
    def __gt__(self, gene):
        if self.chromosome!=gene.chromosome:
            return self.chrom_number>gene.chrom_number
        return self.stop>gene.stop
    def __ge__(self, gene):
        if self.chromosome!=gene.chromosome:
            return self.chrom_number>=gene.chrom_number
        return self.stop>=gene.stop
    def __eq__(self, gene):
        if self.chromosome!=gene.chromosome:
            return False
        return self.stop==gene.stop and self.start==gene.start
    def __neq__(self, gene):
        return not self.__eq__(gene)
    def __len__(self):
        return self.stop-self.start
    def __getitem__(self, item):
        return self.elements[item]
    def get(self, item, default=None):
        if item not in self.elements:
            return default
        return self.elements[item]
    def coding_region_length(self):
        cds_length = 0
        exon_length = 0
        if 'cds' in self.elements:
            for cds in self.elements['cds']:
                cds_length += cds[2]-cds[1]
        if 'exon' in self.elements:
            for exon in self.elements['exon']:
                exon_length += exon[2]-exon[1]
        return max(cds_length, exon_length)
    def __str__(self):
        s = "GENE TRANSCRIPT\n"
        s += "\tChromosome: {}, Start: {},".format(self.chromosome,
                                                   self.start)
        s += " Stop: {}, Direction: {}\n".format(self.stop,
                                                 self.direction)
        if 'cds' in self.elements:
            s += "\tNumber of cds: {}\n".format(len(self.elements['cds']))
        if 'exon' in self.elements:
            s += "\tNumber of exons: {}\n".format(len(self.elements['exon']))
        return s
        
class GFF:
    def __init__(self, filename, queue_size=1):
        if ".gff" in filename: self._filetype = "GFF"
        if ".gtf" in filename: self._filetype = "GTF"
        self._filename = filename
        self._file = open(filename)
        self._Queue = Queue(queue_size)
        self._empty = False
    def get_gene(self):
        if self._empty: return None
        in_gene = False
        while True:
            self._last_line = self._file.tell()
            line = self._file.readline()
            if line=='':
                self._empty = True
                break
            if line[0]=="#": continue
            elements = line.split("\t")
            chrom, _, element_type, start, stop = elements[:5] 
            if in_gene:
                if element_type=="gene" or element_type=="transcript":
                    self._file.seek(self._last_line)
                    break
                gene.add_element(elements)
            else:
                if element_type=="gene" or element_type=="transcript":
                    in_gene = True
                    gene = Gene(elements)
        self._Queue.add(gene)
        return gene
    def get_last_gene(self):
        return self._Queue.get(1)
    def __str__(self):
        return "<--{} file with filename {}-->".format(self._filetype, self._filename)
def test_overlaps(genes):
    overlaps = []
    for i, gene_1 in enumerate(genes):
        for j, gene_2 in enumerate(genes[i:]):
            if i==(j+i): continue
            overlaps.append(gene_1.overlaps(gene_2))
    return overlaps
def get_furthest_gene(genes, empty):
    smaller_than_counts = []
    for i, gene_1 in enumerate(genes):
        counter = 0
        for j, gene_2 in enumerate(genes):
            if i==j: continue
            if empty[i]:
                counter -= 1
                continue
            if gene_1<gene_2: counter += 1
        smaller_than_counts.append(counter)
    return max(enumerate(smaller_than_counts), key=lambda x: x[1])[0]
def exon_overlap(exon1list, exon2list):
    if exon1list==[] or exon2list==[]:
        return 0
    
    n, m = len(exon1list), len(exon2list)
    
    def overlap(exon1, exon2):
        _, start1, stop1, dirc1 = exon1
        _, start2, stop2, dirc2 = exon2
        dirc = dirc1==dirc2
        left = stop2>start1 and start2<=start1
        right = start2<stop1 and stop2>=stop1
        in1 = start1<=start2 and stop1>=stop2
        in2 = start2<=start1 and stop2>=stop1
        if dirc and (left or right or in1 or in2):
            return 1
        else:
            return 0
    O = list() # overlap matrix
    for i in range(n+1):
        O.append(list())
        for j in range(m+1):
            if i==0 or j==0:
                O[i].append(0)
                continue
            ol = overlap(exon1list[i-1], exon2list[j-1])
            O[i].append(max(O[i-1][j], O[i][j-1], O[i-1][j-1]+ol))
    return (O[-1][-1]*2)/float(n+m)
def exon_confidence(exon1list, exon2list):
    if exon1list==[] or exon2list==[]:
        return 0
    cds1 = sum([exon[2]-exon[1] for exon in exon1list])
    cds2 = sum([exon[2]-exon[1] for exon in exon2list])
    
    n, m = len(exon1list), len(exon2list)
    
    def overlap(exon1, exon2):
        _, start1, stop1, dirc1 = exon1
        _, start2, stop2, dirc2 = exon2
        if dirc1!=dirc2: return 0
        area = min(stop1, stop2)-max(start1, start2)
        if area<0: area = 0
        return area
    O = list() # overlap matrix
    for i in range(n+1):
        O.append(list())
        for j in range(m+1):
            if i==0 or j==0:
                O[i].append(0)
                continue
            ol = overlap(exon1list[i-1], exon2list[j-1])
            O[i].append(max(O[i-1][j], O[i][j-1], O[i-1][j-1]+ol))
            
    return (O[-1][-1]*2)/float(cds1+cds2)
    
def exon_overlap_score(genes):
    overlaps_scores = []
    for i, gene_1 in enumerate(genes):
        for j, gene_2 in enumerate(genes[i:]):
            if i==(j+i): continue
            largest_type_1 = max((len(gene_1.get('cds', [])), 'cds'),
                                 (len(gene_1.get('exon', [])), 'exon'),
                                 key=lambda x: x[0])[1]
            largest_type_2 = max((len(gene_2.get('cds', [])), 'cds'),
                                 (len(gene_2.get('exon', [])), 'exon'),
                                 key=lambda x: x[0])[1]
            score = exon_overlap(gene_1.get(largest_type_1, []),
                                 gene_2.get(largest_type_2, []))
            confidence = exon_confidence(gene_1.get(largest_type_1, []),
                                         gene_2.get(largest_type_2, []))
            overlaps_scores.append((score+confidence)/2)
            #overlaps_scores.append(confidence)
            
    return(sum(overlaps_scores)/len(genes), overlaps_scores)
def which_best_overlap(overlaps):
    scores = [0 for _ in range(len(overlaps))]
    c = 0
    for i in range(len(filenames)):
        for j in range(i, len(filenames)):
            if i==j: continue
            scores[i] += overlaps[c]
            scores[j] += overlaps[c]
            c += 1
    return max(enumerate(scores), key=lambda x: x[1])
class GFFwriter:
    def __init__(self, outfile, version="0.0"):
        self._file = open(outfile, "w")
        self._current_version = version
        self._write_header()
        self._transcipt_counter = 1
        self.IDs = {}
        
    def write_gene(self, gene, confidence_level=0, origin=""):
        lines = gene._lines
        #if "ID" in lines[0][-1]:
        #    ID = {key:value for key, value in [l.split("=") for l in lines[0][-1].split(";")]}['ID']
        #else:
        #    ID = lines[0][-1]
        #if ID not in self.IDs:
            #self.IDs[ID] = 0
            
        self._file.write("## Transcript number {}\n".format(self._transcipt_counter))
        self._file.write("## Confidence level for transcript: {}\n".format(confidence_level))
        self._file.write("## transcript origin: {}\n".format(origin))
        self._transcipt_counter += 1
        for line in lines:
            pasteline = "\t".join(line)
            pasteline = pasteline.replace("\n", "")
            pasteline = pasteline.replace("\r", "")
            pasteline = pasteline+"; confidence_level={}\n".format(confidence_level)
            self._file.write(pasteline)
    def _write_header(self):
        self._file.write("## ANNOTATIONcleaner (v{})\n".format(self._current_version))
def Scan_Queues():
    pass
class CandidateBlock:
    pass
def pretty_print(message, colour="reset"):
    base_color = "\033[39m"
    my_color = {"reset": "\033[39m", "green": "\033[32m",
                "cyan": "\033[96m", "blue": "\033[34m",
                "red": "\033[31m", "lightblue": "\033[38;5;74m",
                "orange": "\033[38;5;202m"}
    print(my_color.get(colour, "\033[39m")+message+base_color)
    
if __name__=="__main__":
    _current_version = "0.4.4"
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<gff files or gtf files>\033[39m'''
    parser = OptionParser(usage)
    parser.add_option('-R', type="string", nargs=1, dest="RNAdepth",
                      help="RNA read depth file. <Created by samtools depth -aa reads.bam>")
    parser.add_option('-Q', type="int", nargs=1, dest="QueueSize", default=15,
                      help="Maximum Queue size. Ability to look back to previous transcripts. Larger Q increases runtime and space consumption. Also allows for better comparisons of iso-forms. Default is set to 15.")
    parser.add_option('-c', type="string", nargs=1, dest="Cutoff", default="0.5",
                      help="Confidence level cutoff. Default is set to 0.5.")
    parser.add_option('-s', type="string", nargs=1, dest="SizeLimit", default="30000",
                      help="Gene size limit. Default is set to 30 kb (30000).")
    parser.add_option('-o', type="string", nargs=1, dest="Output", default="output.gtf",
                      help="Output file name. Default: output.gtf")
    pretty_print("======= ANNOTATION MERGER (v{}) =======".format(_current_version), "orange")
    
    options, args = parser.parse_args()
    
    _queue_size = options.QueueSize
    _cut_off = float(options.Cutoff)
    _gene_size_limit = float(options.SizeLimit)
    _outfile_name = options.Output
    pretty_print("Queue size set to {} and Cutoff level set to {}".format(_queue_size, _cut_off), "red")
    pretty_print("Output will be saved into: {}".format(_outfile_name), "red")
    
    pretty_print("OPENING GFF/GTF FILES", "lightblue")
    gene_annotations = []
    filenames = []
    for arg in args:
        if ".gff" in arg or ".gtf" in arg:
            filenames.append(arg)
    
    for filename in filenames:
        pretty_print("\t"+filename)
        gene_annotations.append(GFF(filename))
    pretty_print("PROCESSING GFF/GTF FILES", "lightblue")
    empty = [False, False, False]
    current_genes = [ga.get_gene() for ga in gene_annotations]
    #print(test_overlaps(current_genes))
    transcript_counts = [1 for _ in range(len(gene_annotations))]
    all_overlap_count = 0
    any_overlap_count = 0
    exon_overlaps = []
    passed_genes = 0
    outgff = GFFwriter(_outfile_name)
    display_counter = 1
    display_by = 25000
    pairwise_map = {}
    c = 0
    for i in range(len(filenames)):
        for j in range(i, len(filenames)):
            if i==j: continue
            pairwise_map[c] = (i, j)
            c += 1
            
    ## MAIN SCRIPT
    while (not all(empty)):
        if (display_counter % display_by)==0:
            print("\t{} transcripts read through".format(display_counter))
        overlaps = test_overlaps(current_genes)
        if all(overlaps):
            all_overlap_count += 1
        if any(overlaps):
            any_overlap_count += 1
            score, pair_scores = exon_overlap_score(current_genes)
            if any(score>=_cut_off for score in pair_scores):
                if all(overlaps):
                    bestgene = 0
                    outgff.write_gene(current_genes[bestgene], score, filenames[bestgene])
                else:
                    #which_max, score = max(enumerate(pair_scores), key=lambda x: x[1])
                    #pairs = pairwise_map[which_max]
                    #bestgene = min(pairs)
                    bestgene = which_best_overlap(pair_scores)[0]
                    #bestgene = max(pairs, key=lambda x: len(current_genes[x]))
                    #print pairs
                    #print len(current_genes[pairs[0]]), len(current_genes[pairs[1]])
                    #print bestgene
                    outgff.write_gene(current_genes[bestgene], score, filenames[bestgene])
                    passed_genes += 1
                #if not gene_annotations[bestgene]._empty:
                #    current_genes[bestgene] = gene_annotations[bestgene].get_gene()
                #    while len(current_genes[bestgene])>_gene_size_limit:
                #        current_genes[bestgene] = gene_annotations[bestgene].get_gene()
                #    empty[bestgene] = gene_annotations[bestgene]._empty
                #    transcript_counts[bestgene] += 1
                #    display_counter += 1
                #continue
        lagging_gene = get_furthest_gene(current_genes, empty)
        current_genes[lagging_gene] = gene_annotations[lagging_gene].get_gene()
        while len(current_genes[lagging_gene])>_gene_size_limit:
            current_genes[lagging_gene] = gene_annotations[lagging_gene].get_gene()
        empty[lagging_gene] = gene_annotations[lagging_gene]._empty
        transcript_counts[lagging_gene] += 1
        display_counter += 1
    pretty_print("FILE STATISTICS", "lightblue")
    
    pretty_print("\tAny genes overlap: {}".format(any_overlap_count))
    pretty_print("\tAll genes overlap: {}".format(all_overlap_count))
    pretty_print("\tAverage exon overlap: {}".format(sum(exon_overlaps)/all_overlap_count))
    pretty_print("\tPassed genes: {}".format(passed_genes))
    for filename, tc in zip(filenames, transcript_counts):
         pretty_print("\t{} gene count: {}".format(filename, tc))
    percentage = float(all_overlap_count*len(filenames))/float(sum(transcript_counts))
    pretty_print("\tPercentage all overlap: {} %".format(percentage*100))
    #for ga in gene_annotations:
    #    print(ga._Queue.get(1))
    #for ga in gene_annotations:
    #    print(ga._empty)
    
```

Takes blast results against som database in format: -outfmt '7 qseqid evalue bitscore score stitle' and finds matching queries. It picks the top blast result as a best guess.

``` python
from optparse import OptionParser
import re
class FASTA:
    def __init__(self, filename):
        self._file = open(filename)
        self._empty = False
    def read_sequence(self):
        if self._empty: return None
        inseq = False
        header = ''
        sequence = ''
        while True:
            last_line = self._file.tell()
            line = self._file.readline()
            if line=='':
                self._empty = True
                break
            if line[0]==";": continue
            if line[0]=='>' and inseq:
                self._file.seek(last_line)
                break
            if line[0]=='>' and not inseq:
                header = line[1:-1]
                inseq = True
            if line[0]!='>' and inseq: sequence += line[:-1]
        return(header, sequence)
    def empty(self):
        return self._empty
class BLAST:
    def __init__(self, filename):
        self._file = open(filename)
        self._empty = False
        self._Fields = None
    def read_query(self):
        if self._empty: return None
        in_query = False
        query_name = None
        Matches = []
        while True:
            last_line = self._file.tell()
            line = self._file.readline()
            if line=='':
                self._empty = True
                break
            if line[0]=="#":
                if 'Fields' in line and self._Fields is None:
                    line = line[:-1].split(": ")[-1]
                    self._Fields = line.split(", ")
                continue
            elements = line[:-1].split("\t")
            for i in range(len(elements)):
                element = elements[i]
                try:
                    element = int(element)
                    elements[i] = element
                except:
                    try:
                        element = float(element)
                        elements[i] = element
                    except:
                        pass
            if in_query:
                Matches.append({field:element for field, element in zip(self._Fields, elements)})
                if query_name != elements[0]:
                    self._file.seek(last_line)
                    break
            else:
                query_name = elements[0]
                Matches.append({field:element for field, element in zip(self._Fields, elements)})
                in_query = True
        return query_name, Matches
    def empty(self):
        return self._empty
    
def pretty_print(message, colour="reset"):
    base_color = "\033[39m"
    my_color = {"reset": "\033[39m", "green": "\033[32m",
                "cyan": "\033[96m", "blue": "\033[34m",
                "red": "\033[31m", "lightblue": "\033[38;5;74m",
                "orange": "\033[38;5;202m"}
    print(my_color.get(colour, "\033[39m")+message+base_color)
    
if __name__=="__main__":
    _current_version = "0.1.5"
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<blast result files>\033[39m'''
    parser = OptionParser(usage)
    parser.add_option('-q', type="string", nargs=1, dest="Queries", default=None,
                      help="Input query filename, the sequence used for the blast search. (FASTA format)")
    parser.add_option('-o', type="string", nargs=1, dest="Output", default="output.fasta",
                      help="Output file name. Default: output.fasta")
    parser.add_option('-s', type="string", nargs=1, dest="Summary", default="summary.csv",
                      help="Output summary csv file. Keeps track of the functions given to all transcripts. Default: summary.csv")
    parser.add_option('-e', type="float", nargs=1, dest="Evalue", default=1,
                      help="Evalue cutoff. Default 1")
    pretty_print("======= BLAST ANNOTATER (v{}) =======".format(_current_version), "orange")
    
    options, args = parser.parse_args()
    if options.Queries is None:
        raise "No query file given, please use: -q filename.fasta"
    query = options.Queries
    blast_filename = args[0]
    outfile = options.Output
    summary_file = options.Summary
    eval_cutoff = options.Evalue
    pretty_print("OPENING FILES", "lightblue")
    
    query_file = FASTA(query)
    blast_file = BLAST(blast_filename)
    # araprot regex (\|.+?\|)
    stringmatcher = re.compile(".+?\|.+?\|")
    hexnumbers = re.compile("%..")
    #query1, matches1 = blast.read_query()
    #query2, matches2 = blast.read_query()
    #matches1 = filter(lambda x: x['evalue']<eval_cutoff, matches1)
    #matches2 = filter(lambda x: x['evalue']<eval_cutoff, matches2)
    
    #print len(matches1)
    #print len(matches2)
    #for match in matches1:
    #    print match['subject title']
    #    print stringmatcher.findall(match['subject title'])[0]
    pretty_print("READING BLAST FILE", "lightblue")
    blast_queries = {}
    while not blast_file.empty():
        blast_query, matches = blast_file.read_query()
        blast_queries[blast_query] = matches[0]
    
    pretty_print("READING AND WRITING FASTA FILE", "lightblue")
    
    f = open(outfile, "w")
    s = open(summary_file, "w")
    s.write("transcript;function\n")
    while not query_file.empty():
        gene_name, sequence = query_file.read_sequence()
        gene_name = gene_name.split(" ")[0]
        if gene_name in blast_queries:
            function = stringmatcher.findall(blast_queries[gene_name]['subject title'])[0][:-1]
            if "%" in function:
                hex_values = hexnumbers.findall(function)
                for hex_value in hex_values:
                    function = function.replace(hex_value, hex_value[1:].decode("hex"))
        else:
            function = "NO MATCH"
        s.write(gene_name+";"+function+"\n")
        f.write(">"+gene_name+" | "+function+"\n")
        for i in range((len(sequence)/60)+1):
            f.write(sequence[(i*60):((i+1)*60)])
            f.write("\n")
            
    f.close()
    s.close()
    pretty_print("FINISHED FILES", "lightblue")
    #for match in matches2:
    #    print match
    #print [open(blast_filename).read()]
```

Using the summary csv produced by BLASTannotater.py, it can add function and with query it matched to a gtf file in the attributes column.

``` python
from optparse import OptionParser
def read_summary(summary_file):
    f = open(summary_file)
    summary_data = {}
    f.readline()
    for line in f:
        name, tag = line.split(";", 1)
        summary_data[name] = {}
        tag = tag.replace("\n", "")
        if tag=="NO MATCH":
            summary_data[name]["transcript"] = tag
            summary_data[name]["function"] = "UNKNOWN"
            continue
        transcript, function = tag.split(" | ", 1)
        summary_data[name]["transcript"] = transcript
        summary_data[name]["function"] = function
    return summary_data
def annotate_gff(gff_file, outfile, summary_data):
    gff = open(gff_file)
    gff.readline()
    out = open(outfile, "w")
    for line in gff:
        if line[0]=="#":
            out.write(line)
            continue
        line = line.replace("\n", "")
        line_info = line.split("\t")
        element_type, INFO = line_info[2], line_info[-1]
        if element_type=="transcript":
            INFO = INFO.split(";")
            transcript = INFO[0]
            function = summary_data[transcript]['function']
            matching_query = summary_data[transcript]['transcript']
            compiled_info = {}
            order = []
            for entry in INFO[1:]:
                #print entry
                #print "=" in entry
                if entry == "": continue
                if "=" in entry:
                    key, val = entry.split("=")
                    key.replace(" ", "")
                    compiled_info[key] = val
                    order.append(key)
                else:
                    key, val = entry.split('"', 1)
                    key.replace(" ", "")
                    val = val.replace('"', "")
                    compiled_info[key] = val
                    order.append(key)
            attribute = transcript+"; "
            for key in order:
                attribute += key+' "'+compiled_info[key]+'"; '
            attribute += 'function "'+function+'"; '
            attribute += 'tair_accession "'+matching_query+'"\n'
            
            out.write("\t".join(line_info[:-1])+"\t"+attribute)
        else:
            INFO = INFO.split(";")
            compiled_info = {}
            order = []
            for entry in INFO:
                #print entry
                if entry == "": continue
                if "=" in entry:
                    key, val = entry.split("=")
                    key.replace(" ", "")
                    compiled_info[key] = val
                    order.append(key)
                else:
                    key, val = entry.split('"', 1)
                    key = key.replace(" ", "")
                    val = val.replace('"', "")
                    compiled_info[key] = val
                    order.append(key)
            transcript = compiled_info["transcript_id"]
            function = summary_data[transcript]['function'][:-1]
            matching_query = summary_data[transcript]['transcript']
            attribute = ""
            for key in order:
                attribute += key+' "'+compiled_info[key]+'"; '
            attribute += 'function "'+function+'"; '
            attribute += 'tair_accession "'+matching_query+'"\n'
            
            out.write("\t".join(line_info[:-1])+"\t"+attribute)
            
def pretty_print(message, colour="reset"):
    base_color = "\033[39m"
    my_color = {"reset": "\033[39m", "green": "\033[32m",
                "cyan": "\033[96m", "blue": "\033[34m",
                "red": "\033[31m", "lightblue": "\033[38;5;74m",
                "orange": "\033[38;5;202m"}
    print(my_color.get(colour, "\033[39m")+message+base_color)
if __name__=="__main__":
    _current_version = "0.1"
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<gff file>\033[39m'''
    parser = OptionParser(usage)
    parser.add_option('-o', type="string", nargs=1, dest="Output", default="output.gff", help="Output file name. default: output.gff")
    parser.add_option('-s', type="string", nargs=1, dest="CSV", help="Summary csv file, containing annotations, produced by BLASTannotater.py (Required)")
    pretty_print("======= GFF add info (v{}) =======".format(_current_version), "orange")
    options, args = parser.parse_args()
    if len(args)>0:
        gff_file = args[0]
    else:
        raise "Missing argument of gff_file"
    output = options.Output
    if options.CSV==None:
        raise "Missing summary csv file, please provide it using -s <filename.csv>"
    else:
        summary_file = options.CSV
    pretty_print("READING SUMMARY FILE", "cyan")
    summary_data = read_summary(summary_file)
    pretty_print("ANNOTATING GFF FILE", "cyan")
    annotate_gff(gff_file, output, summary_data)
    
    
```
