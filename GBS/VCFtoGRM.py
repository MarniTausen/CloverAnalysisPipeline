from sys import argv
import numpy as np
import scipy.stats as stats
from optparse import OptionParser

def collectName(string):
    return string.split("/")[-1].split(".")[0]

def getMarker(allele):
    markers = {"0/0": 1.0, "0/1": 0.0, "1/1": -1.0, "./.": None}
    return markers.get(allele.split(":")[0], None)

def readVCFmatrix(filename, n, numpy=True):

    f = open(filename)
    header = ""
    names = []
    matrix = []
    for line in f:
        if line[:2]=="##":
            header += line
            continue
        if line[:1]=="#":
            names = [collectName(name) for name in line.split("\t")[9:]]
            print "Sequence names:", names
            matrix = [[] for i in range(len(names))]
            N = len(names)
            continue
        allele_list = []
        for allele in line.split("\t")[9:]:
            allele_list.append(getMarker(allele))

        if len(filter(lambda x: x!=None, allele_list))>=(N-n):
            for i, allele in enumerate(allele_list):
                if allele!=None:
                    matrix[i].append(allele)
                else:
                    matrix[i].append(0.0)

    if numpy==True:
        matrix = np.array(matrix)

    return matrix, names, header

def replace_column_nans_by_mean(matrix):
    # Set the value of gaps/dashes in each column to be the average of the other values in the column.
    nan_indices = np.where(np.isnan(matrix))
    # Note: bn.nanmean() instead of np.nanmean() because it is a lot(!) faster.
    column_nanmeans = np.nanmean(matrix, axis=0)
    matrix[nan_indices] = np.take(column_nanmeans, nan_indices[1])
    return matrix

def calcGRM(SNP):
    N, M = SNP.shape
    NORM = (SNP-np.mean(SNP, 0))/(np.std(SNP, 0)+0.000000000000001)
    return np.dot(NORM, NORM.T)/M

def help():
    print "====== VCF to GRM (Genetic Relationship Matrix) ====="
    print "Reads in a vcf file as an individual by Marker matrix"
    print "Calculates a GRM using this Marker matrix."
    print "-i <filename.vcf>                 The input vcf file"
    print "-o <filename.csv>                 The output GRM as a csv file"
    print "-n number of missing              How many missing alleles are allowed"

if __name__=="__main__":

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename.csv>")
    parser.add_option('-i', type="string", nargs=1, dest="input", help="<filename.vcf>")
    parser.add_option('-n', type="int", nargs=1, dest="missing", help="number of missing allowed")
    parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
    options, args = parser.parse_args()

    if options.help!=None:
        help()
    if options.input!=None:
        infile = options.input
    else:
        raise "No input file, please include a vcf file -i <filename.vcf>"

    if options.output!=None:
        outfile = options.output
    else:
        outfile = "GRM.csv"

    if options.missing!=None:
        nmissing = int(options.missing)
    else:
        nmissing = 0

    print "READING VCF FILE"
    matrix, names, header = readVCFmatrix(infile, nmissing)

    print "STARTING CALCULATIONS"
    print "Number of individuals", matrix.shape[0]
    M = matrix.shape[1]
    print "Number of Markers", M
    matrix = replace_column_nans_by_mean(matrix)
    GRM = calcGRM(matrix)

    print "SAVING FILE"
    np.savetxt(outfile, GRM, delimiter=", ")
