from sys import argv, stdout

def read_table(filename, sep="\t"):
    f = open(filename)
    header = f.readline().replace("\n", "").split("\t")
    converter = {}
    table = {}
    for i, name in enumerate(header):
        converter[i] = name
        table[name] = []
    for line in f:
        elements = line.replace("\n", "").split(sep)
        for i, element in enumerate(elements):
            table[converter[i]].append(element)
    return table

def write_table(table, columns):
    stdout.write("\t".join(columns)+"\n")
    n = len(table[columns[0]])
    for i in range(n):
        row = []
        for col in columns:
            row.append(str(table[col][i]))
        stdout.write("\t".join(row)+"\n")


def convert_to_popsize(lambda00, mu):
    return (1/float(lambda00))/(2*mu)

def scale_time(time, mu, gentime):
    return (float(time)/mu)*gentime

if __name__=="__main__":
    table = read_table(argv[1])
    mu = float(argv[2])
    gentime = int(argv[3])
    table["scaledPopSize"] = map(lambda x: convert_to_popsize(x, mu), table["lambda_00"])
    table["scaledlefttime"] = map(lambda x: scale_time(x, mu, gentime), table["left_time_boundary"])
    table["scaledrighttime"] = map(lambda x: scale_time(x, mu, gentime), table["right_time_boundary"])

    write_table(table, ["time_index", "scaledlefttime", "scaledrighttime", "scaledPopSize"])
