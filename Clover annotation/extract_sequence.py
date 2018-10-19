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
