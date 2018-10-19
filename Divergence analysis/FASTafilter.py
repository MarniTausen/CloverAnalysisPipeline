import subprocess
from sys import argv, exit

def findlocation(filename, search_term):
    out = subprocess.check_output("grep -b '{}' {}".format(search_term, filename), shell=True)
    matches = out.split("\n")[:-1]
    if len(matches)>1:
        return [int(match.split(":")[0]) for match in matches]
    return int(out.split(":")[0])

def collectsequence(filename, location):
    f = open(filename)
    f.seek(location)
    filebegun = False
    s = ""
    for line in f:
        if filebegun==False:
            if ">" in line:
                filebegun = True
                s += line
        else:
            if ">" in line: break
            s += line
    return s

def filtersequence(filename, location):
    f = open(filename)
    f.seek(location)
    filebegun = False
    for line in f:
        if filebegun==False:
            if ">" in line:
                filebegun = True
                print line,
        else:
            if ">" in line: break
            print line,

if __name__=="__main__":
    search_term, filename = argv[1:3]
    location = findlocation(filename, search_term)
    if type(location)==list:
        seqs = []
        for loc in location:
            #seqs.append([collectsequence(filename, loc), loc])
            filtersequence(filename, loc)
        #location = sorted(seqs, key=lambda x: len(x[0]))[-1][1]
    else:
        filtersequence(filename, location)
