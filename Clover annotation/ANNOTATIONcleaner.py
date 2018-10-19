from ANNOTATIONmerger import *
from math import floor

def overlap(exon1, exon2):
    _, start1, stop1, dirc1 = exon1
    _, start2, stop2, dirc2 = exon2
    dirc = dirc1==dirc2
    left = stop2>start1 and start2<=start1
    right = start2<stop1 and stop2>=stop1
    in1 = start1<=start2 and stop1>=stop2
    in2 = start2<=start1 and stop2>=stop1
    if dirc and (left or right or in1 or in2):
        return 1
    else:
        return 0

def overlap_confidence(exon1, exon2):
    _, start1, stop1, dirc1 = exon1
    _, start2, stop2, dirc2 = exon2
    if dirc1!=dirc2: return 0
    area = min(stop1, stop2)-max(start1, start2)
    if area<0: area = 0
    return area

if __name__=="__main__":
    _current_version = "0.2.7"
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m <main gff file> <gff files to filter by>\033[39m'''

    parser = OptionParser(usage)

    #parser.add_option('-R', type="string", nargs=1, dest="RNAdepth",
    #                  help="RNA read depth file. <Created by samtools depth -aa reads.bam>")
    #parser.add_option('-Q', type="int", nargs=1, dest="QueueSize", default=15,
    #                  help="Maximum Queue size. Ability to look back to previous transcripts. Larger Q increases runtime and space consumption. Also allows for better comparisons of iso-forms. Default is set to 15.")
    parser.add_option('-c', type="string", nargs=1, dest="Cutoff", default="0.5",
                      help="Confidence level cutoff. Default is set to 0.5.")
    parser.add_option('-b', type="string", nargs=1, dest="BinSize", default="25000",
                      help="bin size. Default is set to 25 kb (30000).")
    parser.add_option('-o', type="string", nargs=1, dest="Output", default="output.gtf",
                      help="Output file name. Default: output.gtf")

    pretty_print("======= ANNOTATION CLEANER (v{}) =======".format(_current_version), "orange")
    
    options, args = parser.parse_args()
    
    _cut_off = float(options.Cutoff)
    _outfile_name = options.Output
    _bin_size = int(options.BinSize)
    _queue_size = 1

    pretty_print("Cutoff level set to {} and Bin size set to {}".format(_cut_off, _bin_size), "red")
    pretty_print("Output will be saved into: {}".format(_outfile_name), "red")

    pretty_print("OPENING GFF/GTF FILES", "lightblue")

    gene_annotations = []

    filenames = []
    for arg in args:
        if ".gff" in arg or ".gtf" in arg:
            filenames.append(arg)

    main_gff = GFF(filenames[0])
    pretty_print("\tMain gff file: "+filenames[0])
            
    for filename in filenames[1:]:
        pretty_print("\tAdditional gff files: "+filename)
        gene_annotations.append(GFF(filename))

    pretty_print("READING EXONS FROM THE ADDITIONAL FILES", "lightblue")

    exon_hash = {}
    exon_map = {}

    display_counter = 1
    display_by = 25000

    for ga in gene_annotations:
        while not ga._empty:
            if (display_counter % display_by)==0:
                print("\t{} transcripts read through".format(display_counter))
            gene = ga.get_gene()
            if 'exon' in gene.elements:
                for exon in gene.elements['exon']:
                    string_representation = "%r" % (exon)
                    exon_hash[string_representation] = 0

                    chrom, start, stop, dirc = exon
                    
                    if chrom not in exon_map:
                        exon_map[chrom] = {}
                    if dirc not in exon_map[chrom]:
                        exon_map[chrom][dirc] = {}

                    sbin = int(floor(start/_bin_size))*_bin_size
                    stbin = int(floor(stop/_bin_size))*_bin_size
                    if sbin==stbin:
                        if sbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][sbin] = []
                        exon_map[chrom][dirc][sbin].append(exon)
                    else:
                        if sbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][sbin] = []
                        exon_map[chrom][dirc][sbin].append(exon)
                        if stbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][stbin] = []
                        exon_map[chrom][dirc][stbin].append(exon)
                    
            if 'cds' in gene.elements:
                for exon in gene.elements['cds']:
                    string_representation = "%r" % (exon)
                    exon_hash[string_representation] = 0

                    chrom, start, stop, dirc = exon
                    
                    if chrom not in exon_map:
                        exon_map[chrom] = {}
                    if dirc not in exon_map[chrom]:
                        exon_map[chrom][dirc] = {}

                    sbin = int(floor(start/_bin_size))*_bin_size
                    stbin = int(floor(stop/_bin_size))*_bin_size
                    if sbin==stbin:
                        if sbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][sbin] = []
                        exon_map[chrom][dirc][sbin].append(exon)
                    else:
                        if sbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][sbin] = []
                        exon_map[chrom][dirc][sbin].append(exon)
                        if stbin not in exon_map[chrom][dirc]:
                            exon_map[chrom][dirc][stbin] = []
                        exon_map[chrom][dirc][stbin].append(exon)
            display_counter += 1

    pretty_print("FILTERING MAIN GFF FILE", "lightblue")

    outgff = GFFwriter(_outfile_name, _current_version)
    total_gene_count = 0
    passed_genes = 0
    passed_genes_sub1 = 0
    passed_genes_sub2 = 0

    display_counter = 1
    display_by = 25000

    subgenome1 = set(['chr1', 'chr2', 'chr3', 'chr4',
                      'chr5', 'chr6', 'chr7', 'chr8'])
    subgenome2 = set(['chr9', 'chr10', 'chr11', 'chr12',
                      'chr13', 'chr14', 'chr15', 'chr16'])
    
    while not main_gff._empty:
        if (display_counter % display_by)==0:
                print("\t{} transcripts read through".format(display_counter))
        total_gene_count += 1
        gene = main_gff.get_gene()
        largest_type = max((len(gene.get('cds', [])), 'cds'),
                           (len(gene.get('exon', [])), 'exon'),
                           key=lambda x: x[0])[1]

        exon_count = 0
        exon_overlap_count = 0
        exon_overlap_score = 0
        for exon in gene[largest_type]:
            if "%r" % (exon) in exon_hash:
                exon_count += 1
            chrom, start, stop, dirc = exon
            if chrom in exon_map and dirc in exon_map.get(chrom, {}):
                sbin = int(floor(start/_bin_size))*_bin_size
                stbin = int(floor(stop/_bin_size))*_bin_size
                if sbin==stbin:
                    for test_exon in exon_map[chrom][dirc].get(sbin, []):
                        if overlap(exon, test_exon):
                            exon_overlap_count += 1
                            break
                        #exon_overlap_score += overlap_confidence(exon, test_exon)
                else:
                    for test_exon in exon_map[chrom][dirc].get(sbin, []):
                        if overlap(exon, test_exon):
                            exon_overlap_count += 1
                            break
                    for test_exon in exon_map[chrom][dirc].get(stbin, []):
                        if overlap(exon, test_exon):
                            exon_overlap_count += 1
                            break
        
        confidence = (exon_count+exon_overlap_count)/float(2*len(gene[largest_type]))
        
        if confidence>=_cut_off:
            passed_genes += 1
            if gene.chromosome in subgenome1:
                passed_genes_sub1 += 1
            if gene.chromosome in subgenome2:
                passed_genes_sub2 += 1
            outgff.write_gene(gene, confidence, filenames[0])

        display_counter += 1

    pretty_print("FILE STATISTICS", "lightblue")

    pretty_print("\tTotal gene number: {}".format(total_gene_count))
    pretty_print("\tPassed genes: {}".format(passed_genes))
    pretty_print("\tPassed genes in first subgenome: {}".format(passed_genes_sub1))
    pretty_print("\tPassed genes in second subgenome: {}".format(passed_genes_sub2))
    pretty_print("\tTotal: {}".format(passed_genes_sub1+passed_genes_sub2))
