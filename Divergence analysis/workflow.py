from gwf import Workflow

gwf = Workflow()

#######################################
## Read candidate genes
#######################################

def readgenelist(filename):
    return open(filename).read().split("\n")

def cleangenelist(genelist, title):
    nlist = []
    for gene in genelist:
        if gene=="": continue
        gene = gene.split("To")[1]
        nlist.append(gene)
    return nlist

genelist = readgenelist("Divergentgenes.txt")
#occi_gl = cleangenelist(genelist, "occidentale")
occi_gl = genelist

#############################################
## GMAP database
#############################################

def split_subgenomes(reference, subgenome1, subgenome2):
    """ """
    inputs = [reference]
    outputs = [subgenome1, subgenome2]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    directory = reference.split("/")[0]

    spec = '''
    source activate python2

    python splitreference.py {} {} {} '{}'
    '''.format(reference, subgenome1, subgenome2, ">chr9")

    return inputs, outputs, options, spec

gwf.target_from_template("RepensSplit",
                         split_subgenomes("repens/TrR.v5.fasta",
                                          "repens/TrR.v5.To.fasta",
                                          "repens/TrR.v5.Tp.fasta"))


#############################################
## Making blast databases?
#############################################

#############################################
## genblast
#############################################

def translateSequences(infile, outfile):
    inputs = [infile]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '04:00:00'
    }

    spec = '''
    source activate python2

    python translatefasta.py {} {}
    '''.format(infile, outfile)

    return inputs, outputs, options, spec


def genblast(query_genes, database, outfile, basename=""):
    """ """
    inputs = [query_genes, database]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '04:00:00'
    }

    if basename=="":
        basename = outfile.split("/")[-1].split(".")[0]

    spec = '''
    source activate python2

    cd genBlast

    ./genblast -p genblastg -q ../{query} -t ../{db} -c 1 -r 1 -gff -cdna -o {out}

    python ../cleanfasta.py {out}_1.1c_2.3_s1_0_16_1.DNA > {out}.cleaned.DNA
    python ../FASTafilter.py '\-R1\-' {out}.cleaned.DNA > ../{output}
    '''.format(query=query_genes, db=database, out=basename, output=outfile)

    return inputs, outputs, options, spec


gwf.target_from_template("RepensTogenBlast",
                         genblast("occidentale/divergentgenes.fasta", "repens/TrR.v5.To.fasta", "genblastresults/Togenes.fasta"))
gwf.target_from_template("PalgenBlast",
                         genblast("occidentale/divergentgenes.fasta", "pallescens/final_pallescens.fa", "genblastresults/Palgenes.fasta"))

gwf.target_from_template("Paltranslate",
                         translateSequences("genblastresults/Palgenes.fasta", "genblastresults/Palgenes.prot.fasta"))
gwf.target_from_template("RepensTpgenBlast",
                         genblast("genblastresults/Palgenes.prot.fasta", "repens/TrR.v5.Tp.fasta", "genblastresults/Tpgenes.fasta"))


#############################################
## Generate isolate gene files in occidentale
#############################################

def filter_gene(inputfile, outputfile, searchterm):
    """ """
    inputs = [inputfile]
    outputs = [outputfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '00:20:00'
    }

    spec = '''
    source activate python2

    python FASTafilter.py '{}' {} > {}
    '''.format(searchterm, inputfile, outputfile)

    return inputs, outputs, options, spec

for gene in occi_gl:
    gwf.target_from_template("FilterTo%s" % gene,
                             filter_gene("occidentale/clover.cds.fa",
                                        "genes/To{}.fa".format(gene),
                                         "_"+gene))
    gwf.target_from_template("FilterTp%s" % gene,
                             filter_gene("genblastresults/Palgenes.fasta",
                                        "genes/Tp{}.fa".format(gene),
                                         "_"+gene))
    gwf.target_from_template("FilterTrTo%s" % gene,
                             filter_gene("genblastresults/Togenes.fasta",
                                        "genes/TrTo{}.fa".format(gene),
                                         "_"+gene))
    gwf.target_from_template("FilterTrTp%s" % gene,
                             filter_gene("genblastresults/Tpgenes.fasta",
                                        "genes/TrTp{}.fa".format(gene),
                                         "_"+gene))

#########################################
## Multiple alignment of the genes
#########################################

def join_genes(sequences, outfile, reference_gene):
    inputs = sequences
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '00:20:00'
    }

    spec = '''
    source activate biopython2
    python FASTafilter.py '{}' {} | python joingenes.py'''.format(reference_gene, "occidentale/clover.cds.fa")
    for sequence in sequences:
        spec += " {} ".format(sequence)
    spec += "> {}".format(outfile)

    return inputs, outputs, options, spec

def join_genes_and_translate(sequences, outfile, reference_gene):
    inputs = sequences
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '00:20:00'
    }

    spec = '''
    source activate biopython2
    python FASTafilter.py '{}' {} | python joingenesandtranslate.py'''.format(reference_gene, "occidentale/clover.cds.fa")
    for sequence in sequences:
        spec += " {} ".format(sequence)
    spec += "> {}".format(outfile)

    return inputs, outputs, options, spec


def multiple_align(genes, alignment, settings):
    """ """
    inputs = [genes]
    outputs = [alignment]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }

    spec = '''
    source activate prank

    prank -d={} -o={} {}
    '''.format(genes, ".".join(alignment.split(".")[0:2]), settings)

    return inputs, outputs, options, spec

for gene in occi_gl:
    gwf.target_from_template('makeinput{}'.format(gene),
                             join_genes(['genes/To{}.fa'.format(gene),
                                         'genes/TrTo{}.fa'.format(gene),
                                         'genes/TrTp{}.fa'.format(gene),
                                         'genes/Tp{}.fa'.format(gene)],
                                        "joined_genes/{}.fasta".format(gene), gene))

for gene in occi_gl:
    gwf.target_from_template('alignment{}'.format(gene),
                             multiple_align("joined_genes/{}.fasta".format(gene),
                                            "aligned/{}.best.fas".format(gene), "-codon +F"))


for gene in occi_gl:
    gwf.target_from_template('pepmakeinput{}'.format(gene),
                             join_genes_and_translate(['genes/To{}.fa'.format(gene),
                                         'genes/TrTo{}.fa'.format(gene),
                                         'genes/TrTp{}.fa'.format(gene),
                                         'genes/Tp{}.fa'.format(gene)],
                                        "joined_proteins/{}.fasta".format(gene), gene))

for gene in occi_gl:
    gwf.target_from_template('pepalignment{}'.format(gene),
                             multiple_align("joined_proteins/{}.fasta".format(gene),
                                            "aligned_proteins/{}.best.fas".format(gene), "-protein"))


#########################################
## Calculate dN/dS between the genes
#########################################

def calculate_dNds(alignment, dnds, pdist, sites, dnds_data=""):
    inputs = [alignment]
    outputs = [dnds, pdist, sites, dnds_data]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }

    spec = '''
    source activate python2

    python dNdS.py {} {} {} {} {}
    '''.format(alignment, dnds, pdist, sites, dnds_data)


    return inputs, outputs, options, spec

for gene in occi_gl:
   gwf.target_from_template('dNdS{}'.format(gene),
                            calculate_dNds("aligned/{}.best.fas".format(gene),
                                           "summarydata/{}.dnds.csv".format(gene),
                                           "summarydata/{}.pdist.csv".format(gene),
                                           "summarydata/{}.sites.csv".format(gene),
                                           "summarydata/{}.dnds.info.csv".format(gene)))

def comparePeptides(alignment, pdist, sites):
    inputs = [alignment]
    outputs = [pdist, sites]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }

    spec = '''
    source activate python2

    python comparePeptides.py {} {} {}
    '''.format(alignment, pdist, sites)

    return inputs, outputs, options, spec

for gene in occi_gl:
   gwf.target_from_template('cmpPep{}'.format(gene),
                            comparePeptides("aligned_proteins/{}.best.fas".format(gene),
                                           "summarydata_protein/{}.pdist.csv".format(gene),
                                           "summarydata_protein/{}.sites.csv".format(gene)))

def summary_files(inputfiles, output):
    inputs = inputfiles
    outputs = [output]
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }

    inputstring = ",".join(inputfiles)

    spec = '''
    source activate python2

    python Joinsummary.py "{}" {}
    '''.format(inputstring, output)

    return inputs, outputs, options, spec

blacklist = ["840g41210.1", "0g01140.1", "4488g00060.1", "1089g30080.1", ""]

dNdSfiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    dNdSfiles.append("summarydata/{}.dnds.csv".format(gene))

pdistfiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    pdistfiles.append("summarydata/{}.pdist.csv".format(gene))

sitefiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    sitefiles.append("summarydata/{}.sites.csv".format(gene))

dNdSinfofiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    dNdSinfofiles.append("summarydata/{}.dnds.info.csv".format(gene))

gwf.target_from_template("dNdSsummarybetween",
                         summary_files(dNdSfiles, "summaryfiles/divergent_genes.dnds.csv"))
gwf.target_from_template("pdistsummarybetween",
                         summary_files(pdistfiles, "summaryfiles/divergent_genes.pdist.csv"))
gwf.target_from_template("sitesummarybetween",
                         summary_files(sitefiles, "summaryfiles/divergent_genes.sites.csv"))
gwf.target_from_template("dNdSinfosummarybetween",
                         summary_files(dNdSinfofiles, "summaryfiles/divergent_genes.dnds.info.csv"))


pdistfiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    pdistfiles.append("summarydata_protein/{}.pdist.csv".format(gene))

sitefiles = []
for gene in occi_gl:
    if gene in blacklist: continue
    sitefiles.append("summarydata_protein/{}.sites.csv".format(gene))

gwf.target_from_template("PEPpdistsummarybetween",
                         summary_files(pdistfiles, "summaryfiles/divergent_proteins.pdist.csv"))
gwf.target_from_template("PEPsitesummarybetween",
                         summary_files(sitefiles, "summaryfiles/divergent_proteins.sites.csv"))




############################################
## Sample random genes and run full pipeline
############################################


## LOAD random genes file: random_genes.txt
## Was created by:
## grep '>' occidentale/clover.cds.fa | python createRandomList.py random_genes.txt 30

# randomgenelist = readgenelist("random_genes.txt")

def SplitSampleFasta(filename, outdir, parts):
    inputs = [filename]
    outputs = []
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }

    fdir = "/".join(filename.split("/")[0:-1])

    for i in range(1, parts+1):
        if fdir=="":
            outputs.append(outdir+"/"+filename+".p"+str(i))
        else:
            outputs.append(fdir+"/"+outdir+"/"+filename.split("/")[-1]+".p"+str(i))

    spec = '''
    python splitfasta.py {} {} {}
    '''.format(filename, outdir, parts)

    return inputs, outputs, options, spec

parts = 100

gwf.target_from_template("SplitSamplegenes",
                         SplitSampleFasta("occidentale/clover.protein.fa", "parts", parts))

for i in range(1, parts+1):
    gwf.target_from_template("RepensTogenBlastPart"+str(i),
                         genblast("occidentale/parts/clover.protein.fa.p"+str(i), "repens/TrR.v5.To.fasta",
                                  "genblastresults/parts/Togenes.fasta.p"+str(i), "../occidentale/parts/TrTo.p"+str(i)))
    gwf.target_from_template("PalgenBlastPart"+str(i),
                             genblast("occidentale/parts/clover.protein.fa.p"+str(i), "pallescens/final_pallescens.fa",
                                      "genblastresults/parts/Palgenes.fasta.p"+str(i), "../occidentale/parts/Tp.p"+str(i)))

    gwf.target_from_template("PaltranslatePart"+str(i),
                             translateSequences("genblastresults/parts/Palgenes.fasta.p"+str(i), "genblastresults/parts/Palgenes.prot.fasta.p"+str(i)))
    gwf.target_from_template("RepensTpgenBlastPart"+str(i),
                             genblast("genblastresults/parts/Palgenes.prot.fasta.p"+str(i), "repens/TrR.v5.Tp.fasta",
                                      "genblastresults/parts/Tpgenes.fasta.p"+str(i), "../genblastresults/parts/TrTp.p"))

def join_files(sequences, outfile):
    inputs = sequences
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '00:20:00'
    }

    spec = "cat"
    for sequence in sequences:
        spec += " {} ".format(sequence)
    spec += "> {}".format(outfile)

    return inputs, outputs, options, spec

TrTogenes = []
for i in range(1, parts+1):
    TrTogenes.append("genblastresults/parts/Togenes.fasta.p"+str(i))

Palgenes = []
for i in range(1, parts+1):
    Palgenes.append("genblastresults/parts/Palgenes.fasta.p"+str(i))

TrTpgenes = []
for i in range(1, parts+1):
    TrTpgenes.append("genblastresults/parts/Tpgenes.fasta.p"+str(i))

gwf.target_from_template("AllRepensTo",
                         join_files(TrTogenes, "genblastresults/Togenes.all.fasta"))
gwf.target_from_template("AllRepensPal",
                         join_files(Palgenes, "genblastresults/Palgenes.all.fasta"))
gwf.target_from_template("AllRepensTp",
                         join_files(TrTpgenes, "genblastresults/Tpgenes.all.fasta"))


def efficient_workflow(Togenes, TrTogenes, Tpgenes, Palgenes, gene):
    inputs = [Togenes, TrTogenes, Tpgenes, Palgenes]
    outputs = []
    options = {
        'cores': 2,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '02:00:00'
    }

    spec = '''
    source activate python2

    echo "Filtering the gene from all of the respective species"

    python FASTafilter.py '_{gene}' {Togenes} > /scratch/$SLURM_JOBID/To00.gene.fasta
    python FASTafilter.py '_{gene}' {TrTogenes} > /scratch/$SLURM_JOBID/TrTo.gene.fasta
    python FASTafilter.py '_{gene}' {Tpgenes} > /scratch/$SLURM_JOBID/TrTp.gene.fasta
    python FASTafilter.py '_{gene}' {Palgenes} > /scratch/$SLURM_JOBID/Tp00.gene.fasta

    echo "Joining the gene files"

    python joingenes.py /scratch/$SLURM_JOBID/To00.gene.fasta /scratch/$SLURM_JOBID/TrTo.gene.fasta /scratch/$SLURM_JOBID/TrTp.gene.fasta /scratch/$SLURM_JOBID/Tp00.gene.fasta > /scratch/$SLURM_JOBID/joinedgenes.fasta

    cat /scratch/$SLURM_JOBID/joinedgenes.fasta

    ## Join the genes and translate
    python joingenesandtranslate.py  /scratch/$SLURM_JOBID/To00.gene.fasta /scratch/$SLURM_JOBID/TrTo.gene.fasta /scratch/$SLURM_JOBID/TrTp.gene.fasta /scratch/$SLURM_JOBID/Tp00.gene.fasta > /scratch/$SLURM_JOBID/joinedgenes.prot.fasta

    cat /scratch/$SLURM_JOBID/joinedgenes.prot.fasta

    ## Change environment
    source activate prank

    ## Align genes

    echo "Aligning the genes"

    prank -d=/scratch/$SLURM_JOBID/joinedgenes.fasta -o=/scratch/$SLURM_JOBID/alignedgenes -codon +F

    prank -d=/scratch/$SLURM_JOBID/joinedgenes.prot.fasta -o=/scratch/$SLURM_JOBID/alignedgenes.prot -protein

    cat /scratch/$SLURM_JOBID/alignedgenes.best.fas

    cat /scratch/$SLURM_JOBID/alignedgenes.prot.best.fas

    source activate python2

    ## Calculate measure on the alignments

    python dNdS.py /scratch/$SLURM_JOBID/alignedgenes.best.fas /scratch/$SLURM_JOBID/dnds.csv /scratch/$SLURM_JOBID/pdist.csv /scratch/$SLURM_JOBID/sites.csv /scratch/$SLURM_JOBID/dnds.sites.csv

    python comparePeptides.py /scratch/$SLURM_JOBID/alignedgenes.prot.best.fas /scratch/$SLURM_JOBID/pdist.prot.csv /scratch/$SLURM_JOBID/sites.prot.csv

    cat /scratch/$SLURM_JOBID/dnds.csv
    cat /scratch/$SLURM_JOBID/pdist.csv
    cat /scratch/$SLURM_JOBID/sites.csv
    cat /scratch/$SLURM_JOBID/dnds.sites.csv

    cat /scratch/$SLURM_JOBID/pdist.prot.csv
    cat /scratch/$SLURM_JOBID/sites.prot.csv

    python printcsv.py /scratch/$SLURM_JOBID/dnds.csv "{gene}" all_summaryfiles/all.genes.dnds.csv
    python printcsv.py /scratch/$SLURM_JOBID/pdist.csv "{gene}" all_summaryfiles/all.genes.pdist.csv
    python printcsv.py /scratch/$SLURM_JOBID/sites.csv "{gene}" all_summaryfiles/all.genes.sites.csv
    python printcsv.py /scratch/$SLURM_JOBID/dnds.sites.csv "{gene}" all_summaryfiles/all.genes.dnds.sites.csv

    python printcsv.py /scratch/$SLURM_JOBID/pdist.prot.csv "{gene}" all_summaryfiles/all.genes.pdist.prot.csv
    python printcsv.py /scratch/$SLURM_JOBID/sites.prot.csv "{gene}" all_summaryfiles/all.genes.sites.prot.csv

    '''.format(Togenes=Togenes, TrTogenes=TrTogenes, Tpgenes=Tpgenes, Palgenes=Palgenes, gene=gene)

    return inputs, outputs, options, spec

fullgenelist = readgenelist("allgenes.txt")
#fullgenelist = occi_gl
#fullgenelist = ["56g36650.1", "209g75490.3", "172g65010.1", "912g20040.2", "726g61350.1", "558g23290.1"]

for gene in fullgenelist:
    gwf.target_from_template("Gene"+gene,
                             efficient_workflow("occidentale/clover.cds.fa",
                                                "genblastresults/Togenes.all.fasta",
                                                "genblastresults/Tpgenes.all.fasta",
                                                "genblastresults/Palgenes.all.fasta", gene))
