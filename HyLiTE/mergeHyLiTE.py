from os.path import isfile, join
from glob import glob
from sys import argv
from optparse import OptionParser

def get_unique_samples(filelist):

    sample_names = set()

    for name in filelist:

        name = name.split("Ch.")[-1]
        name = name.split(".")[0]
        sample_names.add(name)

    return sample_names

def which_sample(name, samples):
    for sample in samples:
        if sample+"." in name:
            return sample
    return None

def group_by_sample(filelist, samples):

    new_group = {}
    for sample in samples:
        new_group[sample] = []

    for name in filelist:
        print(name)
        new_group[which_sample(name, samples)].append(name)

    return new_group
        

def joinfiles(filelist, filename):
    newf = open(filename, "w")
    header = None

    for name in filelist:
        f = open(name)
        if header is None:
            header = f.readline()
            newf.write(header)
        else:
            blank = f.readline()
            if header != blank:
                print("Mismatched header! file: {}".format(name))
        newf.write(f.read())
        f.close()

    newf.close()


def process_table(filename):
    f = open(filename)
    header = f.readline()
    header = header.replace("\n", "")
    header = header.split("\t")

    table = {}
    order = []
    
    for line in f:
        line = line.replace("\n", "")
        line = line.split("\t")
        gene = line[0]
        order.append(gene)
        table[gene] = {}
        for key, value in zip(header, line):
            try:
                value = int(value)
            except:
                value = value
            table[gene][key] = value

    return table, order, header

def sum_rows(row1, row2):
    keys = row1.keys()
    newrow = {}
    for key in keys:
        if key=='GENE':
            newrow[key] = row1[key]
            continue
        newrow[key] = row1[key]+row2[key]
    return newrow

def row_sum(row):
    total = 0
    for value in row.values():
        if type(value)==int:
            total += value
    return total
        

def mergefiles(filelist, newfile):
    newf = open(newfile, "w")

    basename = filelist[0]
    table, order, header = process_table(basename)

    for name in filelist[1:]:
        new_table, _, _ = process_table(name)
        for gene in order:
            table[gene] = sum_rows(table[gene], new_table[gene])

    newf.write("\t".join(header)+"\n")

    for gene in order:
        row = []
        for name in header:
            row.append(str(table[gene][name]))
        newf.write("\t".join(row)+"\n")

    newf.close()

def read_protocol_file(filename):
    f = open(filename)
    info = {}
    for line in f:
        line = line.split("\t")
        name, _, sample, _, which = line
        if name not in info:
            info[name] = {}
        info[name][sample] = which
    return info

def extend_expression_information(filename, sample_info, proto_info, outfile):
    expression_table, gene_order, old_header = process_table(filename)

    P1_name = proto_info["P1"]["sample1"].split("/")[-1].split(".")[0]
    P2_name = proto_info["P2"]["sample1"].split("/")[-1].split(".")[0]
    
    new_header = []
    translator = {}
    sample_names = {}
    sample_set = set()
    ch_list = []
    for name in old_header:
        if name=="GENE":
            new_header.append(name)
            translator[name] = name
            continue
        ind, samp = name.split("%")
        origin = proto_info[ind][samp]
        if ind != "Ch":
            species, tissue = origin.split("/")[-1].split(".")[:2]
            compiled_name = "_".join([species, tissue])
            new_header.append(compiled_name)
            translator[compiled_name] = name 
            sample_names[compiled_name] = samp
        else:
            species, tissue = origin.split("/")[-1].split(".")[:2]

            compiled_name_P1 = "_".join([species, tissue, P1_name])
            compiled_name_P2 = "_".join([species, tissue, P2_name])

            new_header.append(compiled_name_P1)
            new_header.append(compiled_name_P2)

            translator[compiled_name_P1] = compiled_name_P1
            translator[compiled_name_P2] = compiled_name_P2

            ch_list.append([compiled_name_P1, compiled_name_P2])

            sample_names[compiled_name_P1] = samp
            sample_names[compiled_name_P2] = samp

            sample_set.add(samp)


    ## ADD NEW ROWS FROM SAMPLE INFO

    for gene in gene_order:
        for chs in ch_list:
            sample = sample_names[chs[0]]
            expression_table[gene][chs[0]] = sample_info[sample][0][gene]["P1"]+sample_info[sample][0][gene]["P1+N"]
            expression_table[gene][chs[1]] = sample_info[sample][0][gene]["P2"]+sample_info[sample][0][gene]["P2+N"]
    
    newf = open(outfile, "w")
        
    newf.write("\t".join(new_header)+"\n")

    for gene in gene_order:
        row = []
        for name in new_header:
            row.append(str(expression_table[gene][translator[name]]))
        newf.write("\t".join(row)+"\n")

    newf.close()

    


def pretty_print(message, colour="reset"):
    base_color = "\033[39m"
    my_color = {"reset": "\033[39m", "green": "\033[32m",
                "cyan": "\033[96m", "blue": "\033[34m",
                "red": "\033[31m", "lightblue": "\033[38;5;74m",
                "orange": "\033[38;5;202m"}
    print(my_color.get(colour, "\033[39m")+message+base_color)
    
if __name__=="__main__":
    _current_version = "0.1.4"
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<additional inputs>\033[39m'''

    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", default="./", help="Output directory")
    parser.add_option('-d', type="string", nargs=1, dest="directory", default="./", help="Input directory of all subsets.")
    parser.add_option('-p', type="string", nargs=1, dest="protocol", default=None, help="Protocol file. (required)")


    pretty_print("======= mergeHyLiTE script (v{}) =======".format(_current_version), "green")

    options, args = parser.parse_args()
    
    directory = options.directory
    outputdir = options.output
    protocol = options.protocol

    if protocol==None:
        raise "No protocol file given, please provide using: -p filename"

    protocol_info = read_protocol_file(protocol)
    
    potential_files = glob(join(directory, "subset*"))
    
    subset_directories = list(filter(lambda x: not isfile(x), potential_files))

    reads = []
    readssummary = []
    expression = []
    snp = []
    snpsummary = []
    runsummary = []

    for sd in subset_directories:
        reads += glob(join(sd, "*.read.txt"))
        readssummary += glob(join(sd, "*.read.summary.txt"))
        expression += glob(join(sd, "*.expression.txt"))
        snp += glob(join(sd, "*.snp.txt"))
        snpsummary += glob(join(sd, "*.snp.summary.txt"))
        runsummary += glob(join(sd, "*.run.summary.txt"))

    sample_names = get_unique_samples(reads)

    print(sample_names)

    readssummary = group_by_sample(readssummary, sample_names)
    reads = group_by_sample(reads, sample_names)

    #print(readssummary)

    for sample in sample_names:
        print("Processing {}".format(sample))
        mergefiles(readssummary[sample], join(outputdir, sample+".read.summary.txt"))

    mergefiles(expression, join(outputdir, "combined.expression.txt"))

    sample_info = {}
    for sample in sample_names:
        sample_info[sample] = process_table(join(outputdir, sample+".read.summary.txt"))

    extend_expression_information(join(outputdir, "combined.expression.txt"), sample_info,
                                  protocol_info, join(outputdir, "completed.expression.txt"))

    

