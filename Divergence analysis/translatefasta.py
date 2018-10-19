from sys import argv

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

def translate(seq):
    aa = []
    for i in range(0, len(seq), 3):
        aa.append(codon_map.get(seq[i:i+3], 'X'))
    return "".join(aa)

def writeFastaAndTranslate(d, filename):
    f = open(filename, "w")
    for k, v in d.items():
        trans = translate(v.upper())
        if trans[-1]=="*": trans = trans[:-1]
        f.write(">%s\n%s\n" % (k, trans))
    f.close()

if __name__=="__main__":
    fasta = readfasta(argv[1])
    writeFastaAndTranslate(fasta, argv[2])
