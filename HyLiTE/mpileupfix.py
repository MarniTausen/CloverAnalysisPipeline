from sys import argv, stdout

def read_and_fix_mpileup(filename, alphabet):
    f = open(filename)
    for line in f:
        line = line.split("\t")
        if line[0]=="retro": continue
        if line[2] not in alphabet: line[2] = "N"
        stdout.write("\t".join(line))
    f.close()

if __name__=="__main__":
    read_and_fix_mpileup(argv[1], "AGTCN")
