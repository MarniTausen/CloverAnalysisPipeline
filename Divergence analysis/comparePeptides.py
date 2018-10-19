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
