from parseVCF import readVCF
from sys import argv, exit
from optparse import OptionParser

def collectReferenceLengths(header):
    contigs = header.split("contig")[1:]
    contigs = [contig.split(">")[0] for contig in contigs]

    lengths = {}
    
    for contig in contigs:
        ID, length = contig.split(",")
        ID = ID.split("=")[-1]
        length = length.split("=")[-1]
        
        lengths[ID] = int(length)
    
    return lengths

def makeBins(referenceSizes, binsize):

    bins = {}
    
    for chromosome, length in referenceSizes.items():
        bins[chromosome] = [{'sites':0, 'SNP':0, 'DP':0, 'AF':0}
                            for _ in range(length/binsize+1)]

    return bins

def countSNPS(data, bins, binsize, depth):

    for i in range(len(data['CHROM'])):
        chrom, pos, alt, dp, af = data['CHROM'][i], data['POS'][i], data['ALT'][i], data['INFO'][i]['DP'], data['INFO'][i].get('AF', 0)

        bpos = pos/binsize

        if dp<depth: continue

        bins[chrom][bpos]['sites'] += 1
        bins[chrom][bpos]['DP'] += dp
        if alt!=".": 
            bins[chrom][bpos]['SNP'] += 1
            try:
                bins[chrom][bpos]['AF'] += af
            except TypeError:
                print "AF was a string:", af

def SNPdensity(bins):
    SNPdens = {}

    for chromosome, counts in bins.items():
        SNPdens[chromosome] = []
        for count in counts:
            if count['sites']==0:
                SNPdens[chromosome].append('NA')
                continue
            SNPdens[chromosome].append(count['SNP']/float(count['sites']))

    return SNPdens

def printSNPdens(SNPdens, binsize):

    print "%20s %20s %20s" % ('Chromosome', 'Location', 'SNP density')
    
    for chromosome, densities in SNPdens.items():
        i = 0
        for density in densities:
            print "%20s %20s %20s" % (chromosome,
                                       str(i*binsize)+"-"+str((i+1)*binsize),
                                       str(density))
            i += 1

def writeCSV(filename, SNPdens, binsize, bins):

    f = open(filename, 'w')
    
    f.write("%s, %s, %s, %s, %s, %s, %s\n" % ('Chromosome', 'Location', 'SNP density', 'Sites', 'SNPs', 'DP', 'AF'))
    
    for chromosome, densities in SNPdens.items():
        for i, density in enumerate(densities):
            sites = bins[chromosome][i]['sites']
            SNP = bins[chromosome][i]['SNP']
            if sites!=0:
                DP = bins[chromosome][i]['DP']/sites
            else: DP = 0
            if SNP!=0:
                AF = bins[chromosome][i]['AF']/float(SNP)
            else: AF = 0
            f.write("%s, %s, %s, %s, %s, %s, %s\n" % (chromosome,
                                                      str(i*binsize)+"-"+str((i+1)*binsize),
                                                      str(density), str(sites),str(SNP),
                                                      str(DP), str(AF)))

    f.close()

def help():
    print "====== SNP density analysis of a gVCF file ====="
    print "Reads in a gvcf file, strips it down and counts"
    print "the occurrences of SNPS versus the number sites sequenced"
    print "-i <filename.vcf>                 The input vcf file"
    print "-o <filename.csv>                 The output csv file"
    print "-b binsize                        The binsize to be calculated from"
    print "-d read depth                     Filter sites less that the minimum read depth"
    print
    exit()

if __name__=="__main__":

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename.csv>")
    parser.add_option('-i', type="string", nargs=1, dest="input", help="<filename.vcf>")
    parser.add_option('-b', type="int", nargs=1, dest="binsize", help="binsize in kb")
    parser.add_option('-d', type="int", nargs=1, dest="depth", help="read depth filter")
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

    if options.binsize!=None:
        binsize = options.binsize*1000
    else:
        raise "Please provide a binsize: -b binsize"

    if options.depth!=None:
        depth = options.depth
    else:
        print "Warning: Was not provided with a minimum depth filter. Therefore it will not filter any sites"
        print "Please use: -d minimum depth, to filter any sites with a low coverage"
        depth = 0

    print "READING DATA:"
    data, header = readVCF(infile, 9)
    print "RETRIEVING REFERENCE LENGTHS"
    reference = collectReferenceLengths(header)
    print "MAKING BINS"
    bins = makeBins(reference, binsize)
    print "COUNTING SNPS"
    countSNPS(data, bins, binsize, depth)
    print "CALCULATING SNP DENSITIES"
    SNPdens = SNPdensity(bins)
    print "WRITING TO CSV"
    writeCSV(outfile, SNPdens, binsize, bins)
