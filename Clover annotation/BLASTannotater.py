from optparse import OptionParser
import re

class FASTA:
    def __init__(self, filename):
        self._file = open(filename)
        self._empty = False

    def read_sequence(self):
        if self._empty: return None
        inseq = False
        header = ''
        sequence = ''
        while True:
            last_line = self._file.tell()
            line = self._file.readline()
            if line=='':
                self._empty = True
                break
            if line[0]==";": continue
            if line[0]=='>' and inseq:
                self._file.seek(last_line)
                break
            if line[0]=='>' and not inseq:
                header = line[1:-1]
                inseq = True
            if line[0]!='>' and inseq: sequence += line[:-1]
        return(header, sequence)

    def empty(self):
        return self._empty

class BLAST:
    def __init__(self, filename):
        self._file = open(filename)
        self._empty = False
        self._Fields = None

    def read_query(self):
        if self._empty: return None
        in_query = False
        query_name = None
        Matches = []
        while True:
            last_line = self._file.tell()
            line = self._file.readline()
            if line=='':
                self._empty = True
                break
            if line[0]=="#":
                if 'Fields' in line and self._Fields is None:
                    line = line[:-1].split(": ")[-1]
                    self._Fields = line.split(", ")
                continue
            elements = line[:-1].split("\t")
            for i in range(len(elements)):
                element = elements[i]
                try:
                    element = int(element)
                    elements[i] = element
                except:
                    try:
                        element = float(element)
                        elements[i] = element
                    except:
                        pass
            if in_query:
                Matches.append({field:element for field, element in zip(self._Fields, elements)})
                if query_name != elements[0]:
                    self._file.seek(last_line)
                    break
            else:
                query_name = elements[0]
                Matches.append({field:element for field, element in zip(self._Fields, elements)})
                in_query = True
        return query_name, Matches

    def empty(self):
        return self._empty

    
def pretty_print(message, colour="reset"):
    base_color = "\033[39m"
    my_color = {"reset": "\033[39m", "green": "\033[32m",
                "cyan": "\033[96m", "blue": "\033[34m",
                "red": "\033[31m", "lightblue": "\033[38;5;74m",
                "orange": "\033[38;5;202m"}
    print(my_color.get(colour, "\033[39m")+message+base_color)
    
if __name__=="__main__":
    _current_version = "0.1.5"
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<blast result files>\033[39m'''

    parser = OptionParser(usage)

    parser.add_option('-q', type="string", nargs=1, dest="Queries", default=None,
                      help="Input query filename, the sequence used for the blast search. (FASTA format)")
    parser.add_option('-o', type="string", nargs=1, dest="Output", default="output.fasta",
                      help="Output file name. Default: output.fasta")
    parser.add_option('-s', type="string", nargs=1, dest="Summary", default="summary.csv",
                      help="Output summary csv file. Keeps track of the functions given to all transcripts. Default: summary.csv")
    parser.add_option('-e', type="float", nargs=1, dest="Evalue", default=1,
                      help="Evalue cutoff. Default 1")

    pretty_print("======= BLAST ANNOTATER (v{}) =======".format(_current_version), "orange")
    
    options, args = parser.parse_args()

    if options.Queries is None:
        raise "No query file given, please use: -q filename.fasta"
    query = options.Queries
    blast_filename = args[0]
    outfile = options.Output
    summary_file = options.Summary
    eval_cutoff = options.Evalue

    pretty_print("OPENING FILES", "lightblue")
    
    query_file = FASTA(query)
    blast_file = BLAST(blast_filename)

    # araprot regex (\|.+?\|)

    stringmatcher = re.compile(".+?\|.+?\|")
    hexnumbers = re.compile("%..")

    #query1, matches1 = blast.read_query()
    #query2, matches2 = blast.read_query()

    #matches1 = filter(lambda x: x['evalue']<eval_cutoff, matches1)
    #matches2 = filter(lambda x: x['evalue']<eval_cutoff, matches2)
    
    #print len(matches1)
    #print len(matches2)

    #for match in matches1:
    #    print match['subject title']
    #    print stringmatcher.findall(match['subject title'])[0]


    pretty_print("READING BLAST FILE", "lightblue")

    blast_queries = {}
    while not blast_file.empty():
        blast_query, matches = blast_file.read_query()
        blast_queries[blast_query] = matches[0]
    
    pretty_print("READING AND WRITING FASTA FILE", "lightblue")
    
    f = open(outfile, "w")
    s = open(summary_file, "w")
    s.write("transcript;function\n")
    while not query_file.empty():
        gene_name, sequence = query_file.read_sequence()
        gene_name = gene_name.split(" ")[0]
        if gene_name in blast_queries:
            function = stringmatcher.findall(blast_queries[gene_name]['subject title'])[0][:-1]
            if "%" in function:
                hex_values = hexnumbers.findall(function)
                for hex_value in hex_values:
                    function = function.replace(hex_value, hex_value[1:].decode("hex"))
        else:
            function = "NO MATCH"
        s.write(gene_name+";"+function+"\n")
        f.write(">"+gene_name+" | "+function+"\n")
        for i in range((len(sequence)/60)+1):
            f.write(sequence[(i*60):((i+1)*60)])
            f.write("\n")
            
    f.close()
    s.close()

    pretty_print("FINISHED FILES", "lightblue")
    #for match in matches2:
    #    print match

    #print [open(blast_filename).read()]
