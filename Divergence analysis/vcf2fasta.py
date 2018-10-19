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
