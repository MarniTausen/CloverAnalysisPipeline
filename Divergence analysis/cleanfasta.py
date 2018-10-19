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
