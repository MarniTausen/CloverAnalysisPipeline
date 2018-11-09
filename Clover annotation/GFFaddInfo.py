from optparse import OptionParser

def read_summary(summary_file):
    f = open(summary_file)
    summary_data = {}
    f.readline()
    for line in f:
        name, tag = line.split(";", 1)
        summary_data[name] = {}
        tag = tag.replace("\n", "")
        if tag=="NO MATCH":
            summary_data[name]["transcript"] = tag
            summary_data[name]["function"] = "UNKNOWN"
            continue
        transcript, function = tag.split(" | ", 1)
        summary_data[name]["transcript"] = transcript
        summary_data[name]["function"] = function
    return summary_data

def annotate_gff(gff_file, outfile, summary_data):
    gff = open(gff_file)
    gff.readline()
    out = open(outfile, "w")
    for line in gff:
        if line[0]=="#":
            out.write(line)
            continue
        line = line.replace("\n", "")
        line_info = line.split("\t")
        element_type, INFO = line_info[2], line_info[-1]
        if element_type=="transcript":
            INFO = INFO.split(";")
            transcript = INFO[0]
            function = summary_data[transcript]['function']
            matching_query = summary_data[transcript]['transcript']

            compiled_info = {}
            order = []
            for entry in INFO[1:]:
                #print entry
                #print "=" in entry
                if entry == "": continue
                if "=" in entry:
                    key, val = entry.split("=")
                    key.replace(" ", "")
                    compiled_info[key] = val
                    order.append(key)
                else:
                    key, val = entry.split('"', 1)
                    key.replace(" ", "")
                    val = val.replace('"', "")
                    compiled_info[key] = val
                    order.append(key)
            attribute = transcript+"; "
            for key in order:
                attribute += key+' "'+compiled_info[key]+'"; '
            attribute += 'function "'+function+'"; '
            attribute += 'tair_accession "'+matching_query+'"\n'
            
            out.write("\t".join(line_info[:-1])+"\t"+attribute)
        else:
            INFO = INFO.split(";")
            compiled_info = {}
            order = []
            for entry in INFO:
                #print entry
                if entry == "": continue
                if "=" in entry:
                    key, val = entry.split("=")
                    key.replace(" ", "")
                    compiled_info[key] = val
                    order.append(key)
                else:
                    key, val = entry.split('"', 1)
                    key = key.replace(" ", "")
                    val = val.replace('"', "")
                    compiled_info[key] = val
                    order.append(key)
            transcript = compiled_info["transcript_id"]
            function = summary_data[transcript]['function'][:-1]
            matching_query = summary_data[transcript]['transcript']

            attribute = ""
            for key in order:
                attribute += key+' "'+compiled_info[key]+'"; '
            attribute += 'function "'+function+'"; '
            attribute += 'tair_accession "'+matching_query+'"\n'
            
            out.write("\t".join(line_info[:-1])+"\t"+attribute)
            

def pretty_print(message, colour="reset"):
    base_color = "\033[39m"
    my_color = {"reset": "\033[39m", "green": "\033[32m",
                "cyan": "\033[96m", "blue": "\033[34m",
                "red": "\033[31m", "lightblue": "\033[38;5;74m",
                "orange": "\033[38;5;202m"}
    print(my_color.get(colour, "\033[39m")+message+base_color)

if __name__=="__main__":
    _current_version = "0.1"
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<gff file>\033[39m'''

    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="Output", default="output.gff", help="Output file name. default: output.gff")
    parser.add_option('-s', type="string", nargs=1, dest="CSV", help="Summary csv file, containing annotations, produced by BLASTannotater.py (Required)")


    pretty_print("======= GFF add info (v{}) =======".format(_current_version), "orange")

    options, args = parser.parse_args()

    if len(args)>0:
        gff_file = args[0]
    else:
        raise "Missing argument of gff_file"

    output = options.Output
    if options.CSV==None:
        raise "Missing summary csv file, please provide it using -s <filename.csv>"
    else:
        summary_file = options.CSV


    pretty_print("READING SUMMARY FILE", "cyan")
    summary_data = read_summary(summary_file)

    pretty_print("ANNOTATING GFF FILE", "cyan")
    annotate_gff(gff_file, output, summary_data)
    
    

