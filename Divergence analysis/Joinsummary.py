from sys import argv

def read_csv_files(filenames):
    full_data = {}
    full_data["ID"] = []
    for filename in filenames:
        f = open(filename)
        for item in full_data:
            full_data[item].append("-")
        full_data["ID"][-1] = filename.split("/")[-1]
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
        for item in row:
            if item not in full_data:
                full_data[item] = ["-" for i in range(len(full_data["ID"]))]
            full_data[item][-1] = row[item]
    return full_data

def write_out_data(filename, full_data):
    f = open(filename, "w")
    sequence_names = full_data["ID"]
    del full_data["ID"]
    keys = sorted(full_data.keys())
    nrows = len(sequence_names)
    f.write("ID")
    for key in keys:
        f.write(","+key)
    f.write("\n")
    for i in range(nrows):
        f.write(sequence_names[i])
        for key in keys:
            f.write(","+full_data[key][i])
        f.write("\n")
    f.close()

if __name__=="__main__":
    filenames = argv[1].split(",")
    combined_data = read_csv_files(filenames)
    write_out_data(argv[2], combined_data)
