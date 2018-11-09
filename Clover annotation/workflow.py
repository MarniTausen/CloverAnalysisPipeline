from gwf import Workflow

gwf = Workflow()

def extract_sequence(reference, sequence_name, output):
    """

    Extract sequence takes in a <reference>, <sequence_name> and <output>
    The <reference> is a fasta file, which must be indexed. (have a .fai) file
    It looks up in the .fai file to find the <sequence_name>.
    If the name is present it reads the sequences from the <reference>,
    and saves it as <output>

    """
    inputs = [reference]
    outputs = [output]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    directory = "/".join(output.split("/")[:-1])

    spec = '''

    mkdir -p {dirc}

    python extract_sequence.py {ref} {seq} > {out}
    '''.format(ref=reference, seq=sequence_name, dirc=directory, out=output)

    return inputs, outputs, options, spec


def repeatmasking(reference, output_reference, directory):
    """
    """
    inputs = [reference]
    outputs = [output_reference]
    options = {
        'cores': 8,
        'memory': '24g',
        'account': 'NChain',
        'walltime': '08:00:00'
    }

    spec = '''
    ../../ANNOTATION/RepeatMasker/bin/RepeatMasker -xsmall -pa {cores} -dir {dirc} -species arabidopsis {ref}
    '''.format(ref=reference, cores=options['cores'], dirc=directory)

    return inputs, outputs, options, spec

def fasta_index(reference):
    inputs = [reference]
    outputs = [reference+".fai"]
    options = {
        'cores': 1,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = '''
    source /com/extra/samtools/1.6.0/load.sh

    samtools faidx {ref}
    '''.format(ref=reference)

    return inputs, outputs, options, spec

gwf.target_from_template("RepensMasking",
                         repeatmasking("repens/TrR.v5.fasta",
                                       "repens/TrR.v5.fasta.masked",
                                       "repens"))

gwf.target_from_template("OccidentaleMasking",
                         repeatmasking("occidentale/To.fasta",
                                       "occidentale/To.fasta.masked",
                                       "occidentale"))

gwf.target_from_template("PallescensMasking",
                         repeatmasking("pallescens/Tp.fasta",
                                       "pallescens/Tp.fasta.masked",
                                       "pallescens"))


gwf.target_from_template("RepensMaskindex",
                         fasta_index("repens/TrR.v5.fasta.masked"))

gwf.target_from_template("OccidentaleMaskindex",
                         fasta_index("occidentale/To.fasta.masked"))

gwf.target_from_template("PallescensMaskindex",
                         fasta_index("pallescens/Tp.fasta.masked"))

def get_sequence_names(reference):
    f = open(reference+".fai")
    sequence_names = []
    for line in f:
        sequence_names.append(line.split("\t")[0])
    return sequence_names

white_clover_chromosomes = get_sequence_names("repens/TrR.v5.fasta")
occidentale_chromosomes = get_sequence_names("occidentale/To.fasta")
pallescens_chromosomes = get_sequence_names("pallescens/Tp.fasta")

n = 10

print("White clover, chromosomes:", len(white_clover_chromosomes))
print("Occidentale, chromosomes:", len(occidentale_chromosomes))
print("Pallescens, chromosomes:", len(pallescens_chromosomes))

for chromosome in white_clover_chromosomes:
    gwf.target_from_template("Repens"+chromosome+"Extract",
                             extract_sequence("repens/TrR.v5.fasta.masked",
                                              chromosome,
                                              "runs/repens/"+chromosome+"/"+chromosome+".fa"))

#for chromosome in occidentale_chromosomes[:n]:
#     gwf.target_from_template("Occidentale"+chromosome+"Extract",
#                              extract_sequence("occidentale/To.fasta.masked",
#                                               chromosome,
#                                               "runs/occidentale/"+chromosome+"/"+chromosome+".fa"))

#for chromosome in pallescens_chromosomes[:n]:
#     gwf.target_from_template("Pallescens"+chromosome+"Extract",
#                              extract_sequence("pallescens/final_pallescens.fa.masked",
#                                               chromosome,
#                                               "runs/pallescens/"+chromosome+"/"+chromosome+".fa"))

gwf.target("OccidentaleSequence", inputs=["occidentale/To.fasta.masked"],
           outputs=["runs/occidentale/occidentalefull/To.fasta"]) << """
mkdir -p runs/occidentale/occidentalefull
cp occidentale/To.fasta.masked runs/occidentale/occidentalefull/To.fasta
"""

gwf.target("PallesecensSequence", inputs=["pallescens/Tp.fasta.masked"],
           outputs=["runs/pallescens/pallescensfull/Tp.fasta"]) << """
mkdir -p runs/pallescens/pallescensfull
cp pallescens/Tp.fasta.masked runs/pallescens/pallescensfull/Tp.fasta
"""

def filter_bam_file(bamfile, chromosome, outfile):
    """

    filter_bam_file uses samtools to read a <bamfile> and read only
    the reads that are mapped to <chromosome>.
    It saves the filtered reads into <outfile>.

    """
    inputs = [bamfile]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    directory = "/".join(outfile.split("/")[:-1])

    spec = '''
    source /com/extra/samtools/1.6.0/load.sh

    mkdir -p {dirc}

    samtools view -b {infile} {chrom} > {out}
    '''.format(infile=bamfile, chrom=chromosome, out=outfile, dirc=directory)

    return inputs, outputs, options, spec



for chromosome in white_clover_chromosomes:
    gwf.target_from_template("Repens"+chromosome+"BamExtract",
                             filter_bam_file("RNAdata/repens/pooled.reads.bam",
                                             chromosome,
                                             "runs/repens/"+chromosome+"/"+chromosome+".bam"))

# for chromosome in occidentale_chromosomes[:n]:
#      gwf.target_from_template("Occidentale"+chromosome+"BamExtract",
#                               filter_bam_file("RNAdata/occidentale/pooled_reads.bam",
#                                               chromosome,
#                                               "runs/occidentale/"+chromosome+"/"+chromosome+".bam"))

# for chromosome in pallescens_chromosomes[:n]:
#      gwf.target_from_template("Pallescens"+chromosome+"BamExtract",
#                               filter_bam_file("RNAdata/pallescens/pooled_reads.bam",
#                                               chromosome,
#                                               "runs/pallescens/"+chromosome+"/"+chromosome+".bam"))


gwf.target("OccidentaleReads", inputs=["RNAdata/occidentale/pooled_reads.bam"],
           outputs=["runs/occidentale/occidentalefull/pooled_reads.bam"]) << """
mkdir -p runs/occidentale/occidentalefull
cp RNAdata/occidentale/pooled_reads.bam runs/occidentale/occidentalefull/pooled_reads.bam
"""

gwf.target("PallesecensReads", inputs=["RNAdata/pallescens/pooled_reads.bam"],
           outputs=["runs/pallescens/pallescensfull/pooled_reads.bam"]) << """
mkdir -p runs/pallescens/pallescensfull
cp RNAdata/pallescens/pooled_reads.bam runs/pallescens/pallescensfull/pooled_reads.bam
"""

def BRAKER_gene_annotation(reference, bamfile, proteinfile, directory, species, outfile, time='12:00:00', cores=4):
    """

    Gene annoation using BRAKER.

    """
    inputs = [reference, bamfile, proteinfile]
    outputs = [outfile]
    options = {
        'cores': cores,
        'memory': '16g',
        'account': 'NChain',
        'walltime': time
    }

    levels = len(directory.split("/"))+1
    bamfile = bamfile.split("/")[-1]
    reference = reference.split("/")[-1]
    proteinfile = "../"*(levels-1)+proteinfile
    
    spec = '''

    cd {directory}
    
    source activate BRAKER

    export GENEMARK_PATH=/faststorage/project/NChain/WHITE_CLOVER/BRAKER_ANNOTATION/gm_et_linux_64/gmes_petap
    export ALIGNMENT_TOOL_PATH=/faststorage/project/NChain/WHITE_CLOVER/BRAKER_ANNOTATION/gth-1.7.1-Linux_x86_64-64bit/bin
    export AUGUSTUS_BIN_PATH=/home/marnit/anaconda3/envs/BRAKER/bin/
    export AUGUSTUS_CONFIG_PATH=/home/marnit/anaconda3/envs/BRAKER/config/
    export AUGUSTUS_SCRIPTS_PATH=/home/marnit/anaconda3/envs/BRAKER/scripts/
    export PERL5LIB=/home/marnit/anaconda3/envs/BRAKER/lib/site_perl/5.26.2
    export PERL5LIB="/home/marnit/anaconda3/envs/BRAKER/lib/site_perl/5.26.2/x86_64-linux-thread-multi:$PERL5LIB"

    export PATH="/home/marnit/anaconda3/envs/BRAKER/bin/:$PATH"

    {calldir}BRAKER_v2.1.0/braker.pl --genome={ref} --prot_seq={prot} --prg=gth --bam={bam} --cores={cores} --species={species} --softmasking --useexisting
    '''.format(directory=directory, calldir="../"*levels, species=species,
               cores=options['cores'], ref=reference, bam=bamfile, prot=proteinfile)

    return inputs, outputs, options, spec


for chromosome in white_clover_chromosomes:
    gwf.target_from_template("Repens"+chromosome+"Annotation",
                             BRAKER_gene_annotation("runs/repens/"+chromosome+"/"+chromosome+".fa",
                                                    "runs/repens/"+chromosome+"/"+chromosome+".bam",
                                                    "MedicagoProtein/Mt.proteins.fa",
                                                    "runs/repens/"+chromosome,
                                                    "repens"+chromosome,
                                                    "runs/repens/"+chromosome+"/braker/repens"+chromosome+"/augustus.hints.gtf"))

# for chromosome in occidentale_chromosomes[:n]:
#     gwf.target_from_template("Occidentale"+chromosome+"Annotation",
#                              BRAKER_gene_annotation("runs/occidentale/"+chromosome+"/"+chromosome+".fa",
#                                                     "runs/occidentale/"+chromosome+"/"+chromosome+".bam",
#                                                     "MedicagoProtein/Mt.proteins.fa",
#                                                     "runs/occidentale/"+chromosome,
#                                                     "occidentale"+chromosome,
#                                                     "runs/occidentale/"+chromosome+"/braker/occidentale"+chromosome+"/augustus.hints.gtf"))

# for chromosome in pallescens_chromosomes[:n]:
#     gwf.target_from_template("Pallescens"+chromosome+"Annotation",
#                              BRAKER_gene_annotation("runs/pallescens/"+chromosome+"/"+chromosome+".fa",
#                                                     "runs/pallescens/"+chromosome+"/"+chromosome+".bam",
#                                                     "MedicagoProtein/Mt.proteins.fa",
#                                                     "runs/pallescens/"+chromosome,
#                                                     "pallescens"+chromosome,
#                                                     "runs/pallescens/"+chromosome+"/braker/pallescens"+chromosome+"/augustus.hints.gtf"))

gwf.target_from_template("OccidentaleAnnotation",
                         BRAKER_gene_annotation("runs/occidentale/occidentalefull/To.fasta",
                                                "runs/occidentale/occidentalefull/pooled_reads.bam",
                                                "MedicagoProtein/Mt.proteins.fa",
                                                "runs/occidentale/occidentalefull",
                                                "occidentale",
                                                "runs/occidentale/occidentalefull/braker/occidentale/augustus.hints.gtf",
                                                "48:00:00",
                                                8))

gwf.target_from_template("PallescensAnnotation",
                         BRAKER_gene_annotation("runs/pallescens/pallescensfull/Tp.fasta",
                                                "runs/pallescens/pallescensfull/pooled_reads.bam",
                                                "MedicagoProtein/Mt.proteins.fa",
                                                "runs/pallescens/pallescensfull",
                                                "pallescens",
                                                "runs/pallescens/pallescensfull/braker/pallescens/augustus.hints.gtf",
                                                "48:00:00",
                                                8))

def join_hints(hints, output):
    """ """
    inputs = hints
    outputs = [output]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = "cat"
    for hint in hints:
        spec += " {}".format(hint)
    spec += " > {}".format(output)

    return inputs, outputs, options, spec


repens_hints = []
for chromosome in white_clover_chromosomes[1:]:
    repens_hints.append("runs/repens/"+chromosome+"/braker/repens"+chromosome+"/augustus.hints.gtf")

gwf.target_from_template("JoinRepensAnnotation",
                         join_hints(repens_hints,
                                    "runs/repens/repens.annotation.gtf"))

# occidentale_hints = []
# for chromosome in occidentale_chromosomes[:n]:
#      occidentale_hints.append("runs/occidentale/"+chromosome+"/braker/occidentale"+chromosome+"/augustus.hints.gtf")

# gwf.target_from_template("JoinOccidentaleAnnoation",
#                          join_hints(occidentale_hints,
#                                     "runs/occidentale/occidentale.annotation.gtf"))


# pallescens_hints = []
# for chromosome in pallescens_chromosomes[:n]:
#      pallescens_hints.append("runs/pallescens/"+chromosome+"/braker/pallescens"+chromosome+"/augustus.hints.gtf")

# gwf.target_from_template("JoinPallescensAnnoation",
#                           join_hints(pallescens_hints,
#                                      "runs/pallescens/pallescens.annotation.gtf"))

gwf.target("OccidentaleFinal", inputs=["runs/occidentale/occidentalefull/braker/occidentale/augustus.hints.gtf"],
           outputs=["runs/occidentale/occidentale.annotation.gtf"]) << """
cp runs/occidentale/occidentalefull/braker/occidentale/augustus.hints.gtf runs/occidentale/occidentale.annotation.gtf
"""

gwf.target("PallesecensFinal", inputs=["runs/pallescens/pallescensfull/braker/pallescens/augustus.hints.gtf"],
           outputs=["runs/pallescens/pallescens.annotation.gtf"]) << """
cp runs/pallescens/pallescensfull/braker/pallescens/augustus.hints.gtf runs/pallescens/pallescens.annotation.gtf
"""


def GMAP_prepare(reference, dbname):
    """ """
    inputs = [reference]
    outputs = [dbname]
    options = {
        'cores': 1,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    directory = dbname.split("/")[0]

    spec = '''
    source /com/extra/gmap/2017-08-15/load.sh

    mkdir -p {dirc}

    gmap_build -d {dbname} {ref} -D $PWD/{dirc}
    '''.format(dbname=dbname.split("/")[-1], ref=reference, dirc=directory)

    return inputs, outputs, options, spec

gwf.target_from_template("OccidentaleDB",
                         GMAP_prepare("occidentale/To.fasta", "gmap/occidentale.gmap"))
gwf.target_from_template("PallescensDB",
                         GMAP_prepare("pallescens/Tp.fasta", "gmap/pallescens.gmap"))


def GMAP_mapping(dbname, genes, gfffile):
    """ """
    inputs = [dbname, genes]
    outputs = [gfffile]
    options = {
        'cores': 8,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '24:00:00'
    }

    directory = dbname.split("/")[0]

    spec = '''
    source /com/extra/gmap/2017-08-15/load.sh

    gmap -d {} -D $PWD/{} -t 8 -f 2 {} > {}

    '''.format(dbname.split("/")[1], directory, genes, gfffile, gfffile, gfffile)

    return inputs, outputs, options, spec

gwf.target_from_template("OccidentaleMedicago",
                         GMAP_mapping("gmap/occidentale.gmap", "MedicagoProtein/Mt4.0v2_Genes.fasta",
                                      "gmap/occidentale.Mt.gff"))
gwf.target_from_template("PallescensMedicago",
                         GMAP_mapping("gmap/pallescens.gmap", "MedicagoProtein/Mt4.0v2_Genes.fasta",
                                      "gmap/pallescens.Mt.gff"))


def gff3plsorting(gff_file, outgff_file):
    """ """
    inputs = [gff_file]
    outputs = [outgff_file]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = '''
    /home/marnit/bin/gff3sort.pl --chr_order natural {infile} > {outfile}
    '''.format(infile=gff_file, outfile=outgff_file)

    return inputs, outputs, options, spec


gwf.target_from_template("OccidentaleMtSort",
                         gff3plsorting("gmap/occidentale.Mt.gff", "gff_files/occidentale.Mt.gff"))

gwf.target_from_template("PallescensMtSort",
                         gff3plsorting("gmap/pallescens.Mt.gff", "gff_files/pallescens.Mt.gff"))

def gffsorting(gff_file, outgff_file):
    """ """
    inputs = [gff_file]
    outputs = [outgff_file]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = '''
    python GFFsort.py -o  {outfile} {infile}
    '''.format(infile=gff_file, outfile=outgff_file)

    return inputs, outputs, options, spec


gwf.target_from_template("OccidentaleOldSort",
                         gffsorting("gff_files/To.v5.gff", "gff_files/To.v5.sorted.gff"))

gwf.target_from_template("PallescensOldSort",
                         gffsorting("gff_files/Tp.v4.gff", "gff_files/Tp.v4.sorted.gff"))

def copy_files(infiles, outfiles):
    """ """
    inputs = infiles
    outputs = outfiles
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = ''
    for inf, outf in zip(infiles, outfiles):
        spec += "cp {} {}\n".format(inf,outf)

    return inputs, outputs, options, spec
    
gwf.target_from_template("MovingFinalGFFfiles",
                         copy_files(["runs/repens/repens.annotation.gtf",
                                     "runs/occidentale/occidentale.annotation.gtf",
                                     "runs/pallescens/pallescens.annotation.gtf"],
                                    ["gff_files/repens.braker.gtf",
                                     "gff_files/occidentale.braker.gtf",
                                     "gff_files/pallescens.braker.gtf"]))

def ANNOTATIONcleaning(main_gff_file, additional_gff_files, output, confidence_level, binsize=25000):
    """ """
    inputs = [main_gff_file]+additional_gff_files
    outputs = [output]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = '''
    source activate python2
python ANNOTATIONmerger.py -h
python ANNOTATIONcleaner.py -o {outfile} -c {conflvl} -b {binsize} {main_gff}'''.format(outfile=output, conflvl=confidence_level,
                                                                                        binsize=binsize, main_gff=main_gff_file)
    for agf in additional_gff_files:
        spec += " {}".format(agf)

    return inputs, outputs, options, spec

gwf.target_from_template("RepensCleaning",
                         ANNOTATIONcleaning("gff_files/repens.braker.gtf",
                                            ["gff_files/TrR.Mt.gff", "gff_files/TrR.v5.sorted.gff"],
                                            "gff_files/repens.final.gtf",
                                            0.25))

gwf.target_from_template("OccidentaleCleaning",
                         ANNOTATIONcleaning("gff_files/occidentale.braker.gtf",
                                            ["gff_files/To.v5.sorted.gff", "gff_files/occidentale.Mt.gff"],
                                            "gff_files/occidentale.final.gtf",
                                            0.25))

gwf.target_from_template("PallescensCleaning",
                         ANNOTATIONcleaning("gff_files/pallescens.braker.gtf",
                                            ["gff_files/Tp.v4.sorted.gff", "gff_files/pallescens.Mt.gff"],
                                            "gff_files/pallescens.final.gtf",
                                            0.25))

def BUSCO(ingff, reference, proteinfile, dbname, report):
    """ """
    inputs = [ingff, reference]
    outputs = [proteinfile, report]
    options = {
        'cores': 4,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    proteinfilename = proteinfile.split("/")[-1]
    
    spec = '''
    source /com/extra/cufflinks/2.2.1/load.sh
    source activate busco

    gffread -y {protfile} -g {ref} {ingff}
    
    cd busco

    run_BUSCO.py -i {protfilename} -o {name} -l ../../busco/sample_data/embryophyta_odb9/ -m proteins -c 4 > {report}
    '''.format(protfile=proteinfile, protfilename=proteinfilename, ref=reference, ingff=ingff, name=dbname, report=report)

    return inputs, outputs, options, spec


gwf.target_from_template("RepensBUSCO",
                         BUSCO("gff_files/repens.final.gtf",
                               "repens/TrR.v5.fasta",
                               "busco/repens.proteins.fa",
                               "RepensFinal",
                               "RepensFinalReport.out"))

gwf.target_from_template("OccidentaleBUSCO",
                         BUSCO("gff_files/occidentale.final.gtf",
                               "occidentale/To.fasta",
                               "busco/occidentale.proteins.fa",
                               "OccidentaleFinal",
                               "OccidentaleFinalReport.out"))

gwf.target_from_template("PallescensBUSCO",
                         BUSCO("gff_files/pallescens.final.gtf",
                               "pallescens/Tp.fasta",
                               "busco/pallescens.proteins.fa",
                               "PallescensFinal",
                               "PallescensFinalReport.out"))
