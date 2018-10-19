#import sys
from optparse import OptionParser
#from matplotlib import pyplot

class Queue:

    def __init__(self, queue_size):
        self._list = list()
        self._max_size = queue_size
        self.full = False

    def add(self, value):
        if self.full:
            self._list.pop(0)
            self._list.append(value)
        else:
            self._list.append(value)
            if len(self._list)==self._max_size:
                self.full = True

    def get(self, element):
        if element>self._max_size: return None
        if len(self._list)<element: return None
        return self._list[-element]

class Gene:

    def __init__(self, information):
        self.elements = {}
        self.new_gene(information)

    def new_gene(self, information):
        chrom, _, element_type, start, stop, _, direction = information[:7]
        self.elements['transcript'] = [chrom, int(start), int(stop), direction]
        self.direction = direction
        self.start = int(start)
        self.stop = int(stop)
        self.chromosome = chrom
        if "_" in chrom:
            self.chrom_number = int(chrom.split("_")[-1])
        else:
            self.chrom_number = int(chrom.split("chr")[-1])
        self._lines = [information]

    def add_element(self, information):
        chrom, _, element_type, start, stop, _, direction = information[:7]
        element_type = element_type.lower()
        if self.direction==".":
            self.direction = direction
        #print self.__contains__([chrom, start, stop, direction])
        if [chrom, int(start), int(stop), direction] in self:
            self._lines.append(information)
            if element_type.lower() not in self.elements:
                self.elements[element_type.lower()] = list()
            self.elements[element_type.lower()].append([chrom, int(start), int(stop), direction])

    def __contains__(self, gene_info):
        if type(gene_info)==Gene:
            return gene_info.direction==self.direction and gene_info.chromosome==self.chromosome and (self.start>=gene_info.start and self.stop<=gene_info.stop)
        if type(gene_info)==list or type(gene_info)==tuple:
            if len(gene_info)==4:
                c, start, stop, d = gene_info
                return d==self.direction and c==self.chromosome and self.start<=start and self.stop>=stop
            else:
                print "WARNING: If list input must be [chromosome, start, stop, direction]"
                return None

    def overlaps(self, gene):
        chrom = gene.chromosome==self.chromosome
        dirc = gene.direction==self.direction
        left = gene.stop>self.start and gene.start<=self.start
        right = gene.start<self.stop and gene.stop>=self.stop
        return (chrom and dirc) and (left or right or (gene in self) or (self in gene))

    def __lt__(self, gene):
        if self.chromosome!=gene.chromosome:
            return self.chrom_number<gene.chrom_number
        return self.stop<gene.stop
    def __le__(self, gene):
        if self.chromosome!=gene.chromosome:
            return self.chrom_number<=gene.chrom_number
        return self.stop<=gene.stop
    def __gt__(self, gene):
        if self.chromosome!=gene.chromosome:
            return self.chrom_number>gene.chrom_number
        return self.stop>gene.stop
    def __ge__(self, gene):
        if self.chromosome!=gene.chromosome:
            return self.chrom_number>=gene.chrom_number
        return self.stop>=gene.stop
    def __eq__(self, gene):
        if self.chromosome!=gene.chromosome:
            return False
        return self.stop==gene.stop and self.start==gene.start
    def __neq__(self, gene):
        return not self.__eq__(gene)
    def __len__(self):
        return self.stop-self.start

    def __getitem__(self, item):
        return self.elements[item]

    def get(self, item, default=None):
        if item not in self.elements:
            return default
        return self.elements[item]

    def coding_region_length(self):
        cds_length = 0
        exon_length = 0
        if 'cds' in self.elements:
            for cds in self.elements['cds']:
                cds_length += cds[2]-cds[1]
        if 'exon' in self.elements:
            for exon in self.elements['exon']:
                exon_length += exon[2]-exon[1]
        return max(cds_length, exon_length)

    def __str__(self):
        s = "GENE TRANSCRIPT\n"
        s += "\tChromosome: {}, Start: {},".format(self.chromosome,
                                                   self.start)
        s += " Stop: {}, Direction: {}\n".format(self.stop,
                                                 self.direction)
        if 'cds' in self.elements:
            s += "\tNumber of cds: {}\n".format(len(self.elements['cds']))
        if 'exon' in self.elements:
            s += "\tNumber of exons: {}\n".format(len(self.elements['exon']))
        return s
        
class GFF:

    def __init__(self, filename, queue_size=1):

        if ".gff" in filename: self._filetype = "GFF"
        if ".gtf" in filename: self._filetype = "GTF"

        self._filename = filename
        self._file = open(filename)
        self._Queue = Queue(queue_size)
        self._empty = False

    def get_gene(self):
        if self._empty: return None
        in_gene = False
        while True:
            self._last_line = self._file.tell()
            line = self._file.readline()
            if line=='':
                self._empty = True
                break
            if line[0]=="#": continue
            elements = line.split("\t")
            chrom, _, element_type, start, stop = elements[:5] 
            if in_gene:
                if element_type=="gene" or element_type=="transcript":
                    self._file.seek(self._last_line)
                    break
                gene.add_element(elements)
            else:
                if element_type=="gene" or element_type=="transcript":
                    in_gene = True
                    gene = Gene(elements)
        self._Queue.add(gene)
        return gene

    def get_last_gene(self):
        return self._Queue.get(1)

    def __str__(self):
        return "<--{} file with filename {}-->".format(self._filetype, self._filename)

def test_overlaps(genes):
    overlaps = []
    for i, gene_1 in enumerate(genes):
        for j, gene_2 in enumerate(genes[i:]):
            if i==(j+i): continue
            overlaps.append(gene_1.overlaps(gene_2))
    return overlaps

def get_furthest_gene(genes, empty):
    smaller_than_counts = []
    for i, gene_1 in enumerate(genes):
        counter = 0
        for j, gene_2 in enumerate(genes):
            if i==j: continue
            if empty[i]:
                counter -= 1
                continue
            if gene_1<gene_2: counter += 1
        smaller_than_counts.append(counter)
    return max(enumerate(smaller_than_counts), key=lambda x: x[1])[0]

def exon_overlap(exon1list, exon2list):

    if exon1list==[] or exon2list==[]:
        return 0
    
    n, m = len(exon1list), len(exon2list)
    
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

    O = list() # overlap matrix
    for i in range(n+1):
        O.append(list())
        for j in range(m+1):
            if i==0 or j==0:
                O[i].append(0)
                continue
            ol = overlap(exon1list[i-1], exon2list[j-1])
            O[i].append(max(O[i-1][j], O[i][j-1], O[i-1][j-1]+ol))

    return (O[-1][-1]*2)/float(n+m)

def exon_confidence(exon1list, exon2list):

    if exon1list==[] or exon2list==[]:
        return 0

    cds1 = sum([exon[2]-exon[1] for exon in exon1list])
    cds2 = sum([exon[2]-exon[1] for exon in exon2list])
    
    n, m = len(exon1list), len(exon2list)
    
    def overlap(exon1, exon2):
        _, start1, stop1, dirc1 = exon1
        _, start2, stop2, dirc2 = exon2
        if dirc1!=dirc2: return 0
        area = min(stop1, stop2)-max(start1, start2)
        if area<0: area = 0
        return area

    O = list() # overlap matrix
    for i in range(n+1):
        O.append(list())
        for j in range(m+1):
            if i==0 or j==0:
                O[i].append(0)
                continue
            ol = overlap(exon1list[i-1], exon2list[j-1])
            O[i].append(max(O[i-1][j], O[i][j-1], O[i-1][j-1]+ol))
            
    return (O[-1][-1]*2)/float(cds1+cds2)
    
def exon_overlap_score(genes):
    overlaps_scores = []
    for i, gene_1 in enumerate(genes):
        for j, gene_2 in enumerate(genes[i:]):
            if i==(j+i): continue
            largest_type_1 = max((len(gene_1.get('cds', [])), 'cds'),
                                 (len(gene_1.get('exon', [])), 'exon'),
                                 key=lambda x: x[0])[1]
            largest_type_2 = max((len(gene_2.get('cds', [])), 'cds'),
                                 (len(gene_2.get('exon', [])), 'exon'),
                                 key=lambda x: x[0])[1]
            score = exon_overlap(gene_1.get(largest_type_1, []),
                                 gene_2.get(largest_type_2, []))
            confidence = exon_confidence(gene_1.get(largest_type_1, []),
                                         gene_2.get(largest_type_2, []))
            overlaps_scores.append((score+confidence)/2)
            #overlaps_scores.append(confidence)
            
    return(sum(overlaps_scores)/len(genes), overlaps_scores)

def which_best_overlap(overlaps):
    scores = [0 for _ in range(len(overlaps))]
    c = 0
    for i in range(len(filenames)):
        for j in range(i, len(filenames)):
            if i==j: continue
            scores[i] += overlaps[c]
            scores[j] += overlaps[c]
            c += 1
    return max(enumerate(scores), key=lambda x: x[1])

class GFFwriter:

    def __init__(self, outfile, version="0.0"):
        self._file = open(outfile, "w")
        self._current_version = version
        self._write_header()
        self._transcipt_counter = 1
        self.IDs = {}
        
    def write_gene(self, gene, confidence_level=0, origin=""):
        lines = gene._lines
        #if "ID" in lines[0][-1]:
        #    ID = {key:value for key, value in [l.split("=") for l in lines[0][-1].split(";")]}['ID']
        #else:
        #    ID = lines[0][-1]
        #if ID not in self.IDs:

            #self.IDs[ID] = 0
            
        self._file.write("## Transcript number {}\n".format(self._transcipt_counter))
        self._file.write("## Confidence level for transcript: {}\n".format(confidence_level))
        self._file.write("## transcript origin: {}\n".format(origin))
        self._transcipt_counter += 1

        for line in lines:
            pasteline = "\t".join(line)
            pasteline = pasteline.replace("\n", "")
            pasteline = pasteline.replace("\r", "")
            pasteline = pasteline+"; confidence_level={}\n".format(confidence_level)
            self._file.write(pasteline)

    def _write_header(self):
        self._file.write("## ANNOTATIONcleaner (v{})\n".format(self._current_version))


def Scan_Queues():
    pass


class CandidateBlock:

    pass

def pretty_print(message, colour="reset"):
    base_color = "\033[39m"
    my_color = {"reset": "\033[39m", "green": "\033[32m",
                "cyan": "\033[96m", "blue": "\033[34m",
                "red": "\033[31m", "lightblue": "\033[38;5;74m",
                "orange": "\033[38;5;202m"}
    print(my_color.get(colour, "\033[39m")+message+base_color)
    
if __name__=="__main__":
    _current_version = "0.4.4"
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<gff files or gtf files>\033[39m'''

    parser = OptionParser(usage)

    parser.add_option('-R', type="string", nargs=1, dest="RNAdepth",
                      help="RNA read depth file. <Created by samtools depth -aa reads.bam>")
    parser.add_option('-Q', type="int", nargs=1, dest="QueueSize", default=15,
                      help="Maximum Queue size. Ability to look back to previous transcripts. Larger Q increases runtime and space consumption. Also allows for better comparisons of iso-forms. Default is set to 15.")
    parser.add_option('-c', type="string", nargs=1, dest="Cutoff", default="0.5",
                      help="Confidence level cutoff. Default is set to 0.5.")
    parser.add_option('-s', type="string", nargs=1, dest="SizeLimit", default="30000",
                      help="Gene size limit. Default is set to 30 kb (30000).")
    parser.add_option('-o', type="string", nargs=1, dest="Output", default="output.gtf",
                      help="Output file name. Default: output.gtf")

    pretty_print("======= ANNOTATION MERGER (v{}) =======".format(_current_version), "orange")
    
    options, args = parser.parse_args()
    
    _queue_size = options.QueueSize
    _cut_off = float(options.Cutoff)
    _gene_size_limit = float(options.SizeLimit)
    _outfile_name = options.Output

    pretty_print("Queue size set to {} and Cutoff level set to {}".format(_queue_size, _cut_off), "red")
    pretty_print("Output will be saved into: {}".format(_outfile_name), "red")
    
    pretty_print("OPENING GFF/GTF FILES", "lightblue")

    gene_annotations = []

    filenames = []
    for arg in args:
        if ".gff" in arg or ".gtf" in arg:
            filenames.append(arg)
    
    for filename in filenames:
        pretty_print("\t"+filename)
        gene_annotations.append(GFF(filename))


    pretty_print("PROCESSING GFF/GTF FILES", "lightblue")

    empty = [False, False, False]
    current_genes = [ga.get_gene() for ga in gene_annotations]
    #print(test_overlaps(current_genes))

    transcript_counts = [1 for _ in range(len(gene_annotations))]
    all_overlap_count = 0
    any_overlap_count = 0
    exon_overlaps = []
    passed_genes = 0

    outgff = GFFwriter(_outfile_name)

    display_counter = 1
    display_by = 25000

    pairwise_map = {}
    c = 0
    for i in range(len(filenames)):
        for j in range(i, len(filenames)):
            if i==j: continue
            pairwise_map[c] = (i, j)
            c += 1
            
    ## MAIN SCRIPT
    while (not all(empty)):
        if (display_counter % display_by)==0:
            print("\t{} transcripts read through".format(display_counter))
        overlaps = test_overlaps(current_genes)
        if all(overlaps):
            all_overlap_count += 1
        if any(overlaps):
            any_overlap_count += 1
            score, pair_scores = exon_overlap_score(current_genes)
            if any(score>=_cut_off for score in pair_scores):
                if all(overlaps):
                    bestgene = 0
                    outgff.write_gene(current_genes[bestgene], score, filenames[bestgene])
                else:
                    #which_max, score = max(enumerate(pair_scores), key=lambda x: x[1])
                    #pairs = pairwise_map[which_max]
                    #bestgene = min(pairs)
                    bestgene = which_best_overlap(pair_scores)[0]
                    #bestgene = max(pairs, key=lambda x: len(current_genes[x]))
                    #print pairs
                    #print len(current_genes[pairs[0]]), len(current_genes[pairs[1]])
                    #print bestgene
                    outgff.write_gene(current_genes[bestgene], score, filenames[bestgene])
                    passed_genes += 1
                #if not gene_annotations[bestgene]._empty:
                #    current_genes[bestgene] = gene_annotations[bestgene].get_gene()
                #    while len(current_genes[bestgene])>_gene_size_limit:
                #        current_genes[bestgene] = gene_annotations[bestgene].get_gene()
                #    empty[bestgene] = gene_annotations[bestgene]._empty
                #    transcript_counts[bestgene] += 1
                #    display_counter += 1
                #continue
        lagging_gene = get_furthest_gene(current_genes, empty)
        current_genes[lagging_gene] = gene_annotations[lagging_gene].get_gene()
        while len(current_genes[lagging_gene])>_gene_size_limit:
            current_genes[lagging_gene] = gene_annotations[lagging_gene].get_gene()
        empty[lagging_gene] = gene_annotations[lagging_gene]._empty
        transcript_counts[lagging_gene] += 1
        display_counter += 1


    pretty_print("FILE STATISTICS", "lightblue")
    
    pretty_print("\tAny genes overlap: {}".format(any_overlap_count))
    pretty_print("\tAll genes overlap: {}".format(all_overlap_count))
    pretty_print("\tAverage exon overlap: {}".format(sum(exon_overlaps)/all_overlap_count))
    pretty_print("\tPassed genes: {}".format(passed_genes))
    for filename, tc in zip(filenames, transcript_counts):
         pretty_print("\t{} gene count: {}".format(filename, tc))
    percentage = float(all_overlap_count*len(filenames))/float(sum(transcript_counts))
    pretty_print("\tPercentage all overlap: {} %".format(percentage*100))

    #for ga in gene_annotations:
    #    print(ga._Queue.get(1))

    #for ga in gene_annotations:
    #    print(ga._empty)
    
