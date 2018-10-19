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
