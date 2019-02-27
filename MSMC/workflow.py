from gwf import Workflow

gwf = Workflow()

def bamCaller(reference, bamfile, vcf_file, mask, mean_depth, chrom):
    inputs = [reference, bamfile]
    outputs = [vcf_file, mask]
    options = {
        'cores': 1,
        'memory': '16g',
        'account': 'NChain',
        'walltime': '24:00:00'
    }

    spec = '''
    source /com/extra/samtools/1.4.1/load.sh
    source /com/extra/bcftools/1.4.1/load.sh

    samtools mpileup -q 20 -Q 20 -C 50 -u -r {chrom} -f {ref} {bam} | bcftools call -c -V indels |
    ./msmc-tools/bamCaller.py {mean_cov} {mask} | gzip -c > {out}
    '''.format(ref=reference, bam=bamfile, mean_cov=mean_depth, mask=mask, out=vcf_file, chrom=chrom)

    return inputs, outputs, options, spec


reference = "/faststorage/project/NChain/WHITE_CLOVER/BRAKER_ANNOTATION/pipeline/repens/TrR.v5.fasta"

individuals = ["ncl-08", "ncl-09", "ncl-10", "ncl-12"]
bamfiles = {"ncl-08": "../ncl-08_files/ncl-08.bam",
            "ncl-09": "../ncl-09_files/ncl-09.bam",
            "ncl-10": "../ncl-10_files/ncl-10.bam",
            "ncl-12": "../ncl-12_files/ncl-12.bam"}
depths = {"ncl-08": 47, "ncl-09": 45,
          "ncl-10": 44, "ncl-12": 55}

chromosomes = ["chr0", "chr1", "chr2", "chr3", "chr4", "chr5",
               "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
               "chr12", "chr13", "chr14", "chr15", "chr16"]

masks = {}
vcf_files = {}

masks_per_individual = {}
vcf_files_per_individual = {}

for individual in individuals:
    if individual not in masks_per_individual:
        masks_per_individual[individual] = []
        vcf_files_per_individual[individual] = []
    for chrom in chromosomes:
        if chrom not in masks:
            masks[chrom] = []
            vcf_files[chrom] = []
        gwf.target_from_template(chrom+"Bamcalling"+individual.split("-")[-1],
                                 bamCaller(reference, bamfiles[individual],
                                           "vcf_files/"+individual+"_"+chrom+".vcf.gz",
                                           "mask_files/"+individual+"_"+chrom+"_mask.bed.gz",
                                           depths[individual], chrom))
        masks[chrom].append("mask_files/"+individual+"_"+chrom+"_mask.bed.gz")
        vcf_files[chrom].append("vcf_files/"+individual+"_"+chrom+".vcf.gz")
        masks_per_individual[individual].append("mask_files/"+individual+"_"+chrom+"_mask.bed.gz")
        vcf_files_per_individual[individual].append("vcf_files/"+individual+"_"+chrom+".vcf.gz")


def generate_multihetsep(masks, vcf_files, output):
    inputs = masks+vcf_files
    outputs = [output]
    options = {
        'cores': 1,
        'memory': '16g',
        'account': 'NChain',
        'walltime': '24:00:00'
    }

    spec = "./msmc-tools/generate_multihetsep.py "
    for mask in masks:
        spec += "--mask={} ".format(mask)
    for vcf_file in vcf_files:
        spec += "{} ".format(vcf_file)
    spec += "> {}".format(output)

    return inputs, outputs, options, spec

multihetsep_haplo8 = []

for chrom in chromosomes:
    gwf.target_from_template("generate_multihetsep_{}_hap8".format(chrom),
                             generate_multihetsep(masks[chrom], vcf_files[chrom],
                                                  "8haplo/clover_{}.mhs".format(chrom)))
    multihetsep_haplo8.append("8haplo/clover_{}.mhs".format(chrom))

haplo6_ind = ["ncl-09", "ncl-10", "ncl-12"]
mhs_haplo6 = []

temp_masks = {}
temp_vcf = {}

for individual in haplo6_ind:
    for chrom, mask, vcf in zip(chromosomes, masks_per_individual[individual], vcf_files_per_individual[individual]):
        if chrom not in temp_masks:
            temp_masks[chrom] = []
            temp_vcf[chrom] = []
        temp_masks[chrom].append(mask)
        temp_vcf[chrom].append(vcf)

for chrom in chromosomes:
    gwf.target_from_template("generate_multihetsep_{}_hap6".format(chrom),
                             generate_multihetsep(temp_masks[chrom], temp_vcf[chrom],
                                                  "6haplo/clover_{}.mhs".format(chrom)))
    mhs_haplo6.append("6haplo/clover_{}.mhs".format(chrom))


haplo4 = {"1": ["ncl-09", "ncl-10"], "2": ["ncl-10", "ncl-12"],
          "3" : ["ncl-09", "ncl-12"]}
mhs_haplo4 = {}

for case in haplo4:

    temp_masks = {}
    temp_vcf = {}

    mhs_haplo4[case] = []

    for individual in haplo4[case]:
        for chrom, mask, vcf in zip(chromosomes, masks_per_individual[individual], vcf_files_per_individual[individual]):
            if chrom not in temp_masks:
                temp_masks[chrom] = []
                temp_vcf[chrom] = []
            temp_masks[chrom].append(mask)
            temp_vcf[chrom].append(vcf)

    for chrom in chromosomes:
        gwf.target_from_template("generate_multihetsep_{}_hap4_{}".format(chrom, case),
                                 generate_multihetsep(temp_masks[chrom], temp_vcf[chrom],
                                                      "4haplo/clover_{}_{}.mhs".format(chrom, case)))
        mhs_haplo4[case].append("4haplo/clover_{}_{}.mhs".format(chrom, case))

haplo2 = {"1": ["ncl-08"], "2": ["ncl-09"],
          "3" : ["ncl-10"], "4": ["ncl-12"]}
mhs_haplo2 = {}

for case in haplo2:

    temp_masks = {}
    temp_vcf = {}

    mhs_haplo2[case] = []

    for individual in haplo2[case]:
        for chrom, mask, vcf in zip(chromosomes, masks_per_individual[individual], vcf_files_per_individual[individual]):
            if chrom not in temp_masks:
                temp_masks[chrom] = []
                temp_vcf[chrom] = []
            temp_masks[chrom].append(mask)
            temp_vcf[chrom].append(vcf)

    for chrom in chromosomes:
        gwf.target_from_template("generate_multihetsep_{}_hap2_{}".format(chrom, case),
                                 generate_multihetsep(temp_masks[chrom], temp_vcf[chrom],
                                                      "2haplo/clover_{}_{}.mhs".format(chrom, case)))
        mhs_haplo2[case].append("2haplo/clover_{}_{}.mhs".format(chrom, case))


def msmc(multihetsep, output_name):
    inputs = multihetsep
    outputs = [output_name+".log", output_name+".loop.txt",
               output_name+".final.txt"]
    options = {
        'cores': 8,
        'memory': '68g',
        'account': 'NChain',
        'walltime': '24:00:00'
    }

    spec = "./msmc --fixedRecombination -r 0.15 -o {} -t {} ".format(output_name, options['cores'])
    for mhs in multihetsep:
        spec += "{} ".format(mhs)

    return inputs, outputs, options, spec

gwf.target_from_template("msmchaplo8",
                         msmc(multihetsep_haplo8, "clover_4_ind"))

gwf.target_from_template("msmchaplo6",
                         msmc(mhs_haplo6, "clover_3_ind"))

for case in haplo4:
    gwf.target_from_template("msmchaplo4_{}".format(case),
                             msmc(mhs_haplo4[case], "clover_2_ind_{}".format(case)))

for case in haplo2:
    gwf.target_from_template("msmchaplo2_{}".format(case),
                             msmc(mhs_haplo2[case], "clover_1_ind_{}".format(case)))

