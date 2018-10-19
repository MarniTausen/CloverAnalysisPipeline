from sys import argv

def readfasta(filename):
    sequences = open(filename).read().split(">")[1:]
    sequences = [seq.split("\n", 1) for seq in sequences]
    return {seq[0]:seq[1].replace("\n", "") for seq in sequences}

def writeFasta(d, filename):
    f = open(filename, "w")
    for k, v in d.items():
        f.write(">%s\n%s\n" % (k, v))
    f.close()

if __name__=="__main__":
    filename = argv[1]
    fasta = readfasta(argv[1])
    outdir = "/".join(filename.split("/")[0:-1])+"/"+argv[2]
    if "/"==outdir[0]: outdir = outdir[1:]
    splits = int(argv[3])
    gene_names = list(fasta.keys())
    n = len(gene_names)
    size = n//splits
    for i in range(0, splits):
        ndict = {}
        for gn in gene_names[i*size:(i+1)*size]:
            ndict[gn] = fasta[gn]
        writeFasta(ndict, outdir+"/"+filename.split("/")[-1]+".p"+str(i+1))
