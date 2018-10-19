from sys import argv

def make_new_reference_files(filename, sub1, sub2, divider=">chr9"):
    genomes = open(filename).read().split(divider)
    f = open(sub1, "w")
    f.write(genomes[0])
    f.close()
    f = open(sub2, "w")
    f.write(">chr9"+genomes[1])
    f.close()

if __name__=="__main__":
    make_new_reference_files(argv[1], argv[2], argv[3], argv[4])
