from sys import argv
import fcntl

combinations = ['To|Tp', 'To|TrTo', 'To|TrTp', 'Tp|TrTo', 'Tp|TrTp', 'TrTo|TrTp']

def read_csv_and_print(filename, ID, outfile):
    f = open(filename)
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
    f.close()

    s = "C"+ID
    for comb in combinations:
        s += ",{}".format(row.get(comb, "NA"))
    s += "\n"

    g = open(outfile, "a")
    fcntl.flock(g, fcntl.LOCK_EX)
    g.write(s)
    fcntl.flock(g, fcntl.LOCK_UN)
    g.close()

if __name__=="__main__":
    read_csv_and_print(argv[1], argv[2], argv[3])
