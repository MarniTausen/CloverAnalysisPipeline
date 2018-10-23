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

