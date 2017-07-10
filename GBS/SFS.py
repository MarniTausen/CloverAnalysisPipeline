from sys import argv
from parseVCF import readVCF
from optparse import OptionParser

def makeSFS(data, depth, samples):

    sfs = [0 for i in range(0, samples*2)]

    for info in data:
        if info['DP']<depth: continue
        if info['AN']<samples: continue
        try:
            sfs[int(info['AC'])] += 1
        except:
            pass

    return sfs

def writeSFS(SFS, outfile):

    f = open(outfile, "w")
    f.write("AN, Count\n")
    for i, c in enumerate(SFS):
        f.write("%i, %i\n" % (i, c))
    f.close()

def help():
    print "====== S Site Frequency spectrum alysis of vcf file ====="
    print "Reads in a vcf file, strips it down and counts"
    print "the occurrences of SNPS versus the number sites sequenced"
    print "-i <filename.vcf>                 The input vcf file"
    print "-o <filename.csv>                 The output csv file"
    print "-d depth                          Read depth to filter by"
    print "-s samples                        Number of samples"

if __name__=="__main__":

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename.csv>")
    parser.add_option('-i', type="string", nargs=1, dest="input", help="<filename.vcf>")
    parser.add_option('-d', type="int", nargs=1, dest="depth", help="read depth filter")
    parser.add_option('-s', type="int", nargs=1, dest="samples", help="number of samples")
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
        outfile = "output.csv"

    if options.depth!=None:
        depth = options.depth
    else:
        print "Warning: Was not provided with a minimum depth filter. Therefore it will not filter any sites"
        print "Please use: -d minimum depth, to filter any sites with a low coverage"
        depth = 0

    if options.samples!=None:
        samples = options.samples
    else:
        raise "Sample size unknown. Please use: -s samples, to filter sites by the number of samples."

    INFO = readVCF(infile, 9)[0]['INFO']
    sfs = makeSFS(INFO, depth, samples)
    writeSFS(sfs, outfile)



