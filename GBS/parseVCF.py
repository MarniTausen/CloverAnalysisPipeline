from sys import argv

def simplify_name(name):
    if "/" in name:
        return name.split("/")[-1].split(".")[0]
    return name

def parseINFO(info):
    stats = {}
    for statistic in info.split(";")[:6]:
        name, value = statistic.split("=")
        try:
            stats[name] = float(value)
        except ValueError:
            stats[name] = value
    return stats

def readVCF(filename, nc=6):
    ''' 
    Reads in a VCF file while discarding most of the data.
    Only reads in, Chromosome (CHROM), position (POS), reference nucleotide (REF),
    alternative nucleotide (ALT), and mapping quality (QUAL). 
    '''

    f = open(filename)
    header = ""
    names = []
    data = {}
    options = {'CHROM': lambda x: data['CHROM'].append(x),
               'POS': lambda x: data['POS'].append(int(x)),
               'ID': lambda x: 1,
               'REF': lambda x: data['REF'].append(x),
               'ALT': lambda x: data['ALT'].append(x),
               'QUAL': lambda x: data['QUAL'].append(float(x)),
               'INFO': lambda x: data['INFO'].append(parseINFO(x)),
               'FORMAT': lambda x: 1,
               'FILTER': lambda x: 1
    }
    for line in f:
        if line[:2]=="##":
            header += line
            continue
        if line[:1]=="#":
            names = line.split("\t")[:nc]
            names = [simplify_name(name) for name in names]
            names[0] = names[0][1:]
            for name in names:
                data[name] = []
            print names
            continue
        for i, point in enumerate(line.split("\t")[:nc]):
            name = names[i]
            if name not in options:
                data[name].append(point)
            else: 
                options[name](point)

    del data['ID']
    
    return data, header


if __name__=="__main__":
    data, header = readVCF(argv[1], int(argv[2]))
    print data
