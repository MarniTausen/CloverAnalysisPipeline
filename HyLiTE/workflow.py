from gwf import Workflow

gwf = Workflow()

## Reference is a fasta file of the genes.

def bowtie2_index(reference, reference_name):
    inputs = [reference]
    outputs = [reference_name+".1.bt2"]
    options = {
        'cores': 1,
        'memory': '16g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = '''
    source activate HyLiTE

    bowtie2-build {ref} {ref_n}
    '''.format(ref=reference, ref_n=reference_name)

    return inputs, outputs, options, spec

def star_index(reference, output):
    inputs = [reference]
    outputs = [output]
    options = {"cores":8, "memory":"64g", "account":"NChain", "walltime": "12:00:00"}

    directory = "/".join(output.split("/")[:-1])

    spec = """
    source /com/extra/STAR/2.5.2b/load.sh
    STAR --runMode genomeGenerate --runThreadN 8 --genomeDir {dir} --genomeFastaFiles {ref} --limitGenomeGenerateRAM=64000000000 --genomeSAindexNbases 3
    """.format(ref=reference, dir=directory)

    return inputs, outputs, options, spec


reference_file = "./references/To/To.v5.gDNA.fasta"
index_ref_file = "./references/To/To.ref"

raw_gDNA = {"To": ["/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_occidentale_180bp_1.fastq.gz",
                   "/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_occidentale_180bp_2.fastq.gz"],
            "Tp": ["/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_pallescens_180_1.fastq.gz",
                   "/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_pallescens_180_2.fastq.gz"],
            "TrR": ["/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_repens_180bp_1.fastq.gz",
                    "/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_repens_180bp_2.fastq.gz"]}

raw_RNA_old = {"To": {"floral": ["../BRAKER_ANNOTATION/pipeline/RNAdata/occidentale/To_F2_pooled_1.fastq",
                             "../BRAKER_ANNOTATION/pipeline/RNAdata/occidentale/To_F2_pooled_2.fastq"],
                  "leaf": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_1_1.fastq",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_1_2.fastq"],
                  "stolon": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_3_1.fastq",
                             "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_3_2.fastq"],
                  "root": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_6_1.fastq",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_6_2.fastq"]},
           "Tp": {"floral": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_F1_pooled_1.fastq",
                             "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_F1_pooled_2.fastq"],
                  "leaf": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_L4_1_1.fastq",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_L4_1_2.fastq"],
                  "stolon": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/T_pal_stolon_1_1_9_1.fastq",
                             "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/T_pal_stolon_1_1_9_2.fastq"],
                  "root": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_root_3e1_47_1.fastq",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_root_3e1_47_2.fastq"]},
           "TrR": {"floral": ["../BRAKER_ANNOTATION/pipeline/RNAdata/repens/floral/TRfloral_1.fastq",
                              "../BRAKER_ANNOTATION/pipeline/RNAdata/repens/floral/TRfloral_2.fastq"],
                   "leaf": ["../BRAKER_ANNOTATION/pipeline/RNAdata/repens/leaf/TRleaf_1.fastq",
                            "../BRAKER_ANNOTATION/pipeline/RNAdata/repens/leaf/TRleaf_2.fastq"],
                   "stolon": ["../BRAKER_ANNOTATION/pipeline/RNAdata/repens/stolon/TRstolon_1.fastq",
                              "../BRAKER_ANNOTATION/pipeline/RNAdata/repens/stolon/TRstolon_2.fastq"],
                   "root": ["../BRAKER_ANNOTATION/pipeline/RNAdata/repens/root/TRroot_1.fastq",
                            "../BRAKER_ANNOTATION/pipeline/RNAdata/repens/root/TRroot_2.fastq"]}}



raw_RNA = {"To": {"YL1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_YL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_YL_2.fq.gz"],
                  "YL2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_YL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_YL_2.fq.gz"],
                  "YL3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_YL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_YL_2.fq.gz"],
                  "OL1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_OL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_OL_2.fq.gz"],
                  "OL2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_OL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_OL_2.fq.gz"],
                  "OL3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_OL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_OL_2.fq.gz"],
                  "St1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_St_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_St_2.fq.gz"],
                  "St2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_St_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_St_2.fq.gz"],
                  "St3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_St_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_St_2.fq.gz"],
                  "R1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_R_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_R_2.fq.gz"],
                  "R2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_R_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_R_2.fq.gz"],
                  "R3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_R_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_R_2.fq.gz"]},
           "Tp": {"YL1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_YL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_YL_2.fq.gz"],
                  "YL2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_YL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_YL_2.fq.gz"],
                  "YL3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_YL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_YL_2.fq.gz"],
                  "OL1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_OL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_OL_2.fq.gz"],
                  "OL2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_OL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_OL_2.fq.gz"],
                  "OL3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_OL_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_OL_2.fq.gz"],
                  "St1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_St_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_St_2.fq.gz"],
                  "St2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_St_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_St_2.fq.gz"],
                  "St3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_St_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_St_2.fq.gz"],
                  "R1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_R_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_R_2.fq.gz"],
                  "R2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_R_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_R_2.fq.gz"],
                  "R3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_R_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_R_2.fq.gz"]},
           "S9": {"F1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_108_F_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_108_F_2.fq.gz"],
                  "F2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_60_F_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_60_F_2.fq.gz"],
                  "F3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_86_F_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_86_F_2.fq.gz"],
                  "L1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_108_L_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_108_L_2.fq.gz"],
                  "L2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_60_L_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_60_L_2.fq.gz"],
                  "L3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_86_L_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_86_L_2.fq.gz"],
                  "St1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_108_St_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_108_St_2.fq.gz"],
                  "St2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_60_St_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_60_St_2.fq.gz"],
                  "St3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_86_St_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_86_St_2.fq.gz"],
                  "R1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_108_R_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_108_R_2.fq.gz"],
                  "R2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_60_R_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_60_R_2.fq.gz"],
                  "R3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_86_R_1.fq.gz",
                         "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S9_86_R_2.fq.gz"]},
           "S10": {"St1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_1_St_1.fq.gz",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_1_St_2.fq.gz"],
                   "St2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_2_St_1.fq.gz",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_2_St_2.fq.gz"],
                   "St3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_3_St_1.fq.gz",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_3_St_2.fq.gz"],
                   "YL1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_1_YL_1.fq.gz",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_1_YL_2.fq.gz"],
                   "YL2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_2_YL_1.fq.gz",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_2_YL_2.fq.gz"],
                   "YL3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_3_YL_1.fq.gz",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_3_YL_2.fq.gz"],
                   "OL1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_1_OL_1.fq.gz",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_1_OL_2.fq.gz"],
                   "OL2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_2_OL_1.fq.gz",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_2_OL_2.fq.gz"],
                   "OL3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_3_OL_1.fq.gz",
                           "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_3_OL_2.fq.gz"],
                   "R1": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_1_R_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_1_R_2.fq.gz"],
                   "R2": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_2_R_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_2_R_2.fq.gz"],
                   "R3": ["/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_3_R_1.fq.gz",
                          "/home/marnit/NChain/faststorage/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/S10_3_R_2.fq.gz"]}}

all_species_gDNA = ["To", "Tp", "TrR"]
all_species_RNA = ["To", "Tp", "S9", "S10"]
tissues_old = ["floral", "leaf", "stolon", "root"]

tissues = {"To": ["St1", "St2", "St3", "R1", "R2", "R3", "YL1", "YL2", "YL3", "OL1", "OL2", "OL3"],
           "Tp": ["St1", "St2", "St3", "R1", "R2", "R3", "YL1", "YL2", "YL3", "OL1", "OL2", "OL3"],
           "S9": ["F1", "F2", "F3", "St1", "St2", "St3", "R1", "R2", "R3", "L1", "L2", "L3"],
           "S10": ["St1", "St2", "St3", "R1", "R2", "R3", "YL1", "YL2", "YL3", "OL1", "OL2", "OL3"]}


gwf.target_from_template("ToIndex",
                         bowtie2_index(reference_file, index_ref_file))

def star_mapping(read1, output, output_prefix, genomeDir):
    inputs = [read1, genomeDir+"/genomeParameters.txt"]
    outputs = [output]
    options = {
	"cores": 8,
	"memory": "64g",
	"account": "NChain",
	"walltime": "12:00:00"}
    OFilePrefix = output_prefix

    spec = """
    source /com/extra/STAR/2.5.2b/load.sh
    STAR --runThreadN 8 --genomeDir {dir} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {of} --readFilesIn {r1} --limitGenomeGenerateRAM=64000000000""".format(r1 = read1, of=output_prefix, dir=genomeDir)

    if read1.split(".")[-1]=="gz":
        spec += " --readFilesCommand zcat"

    return inputs, outputs, options, spec

def bowtie2_mapping(read1, output, genomeDir):
    inputs = [read1, genomeDir+".1.bt2"]
    outputs = [output]
    options = {
	"cores": 8,
	"memory": "24g",
	"account": "NChain",
	"walltime": "12:00:00"}

    spec = """
    source activate HyLiTE

    bowtie2 -N 1 -x {dir} -U {reads} -S {sam} --local -p {cores}
    """.format(reads = read1, sam=output, dir=genomeDir, cores=options["cores"])

    return inputs, outputs, options, spec


for species in all_species_gDNA:
    for i, readfile in enumerate(raw_gDNA[species]):
        gwf.target_from_template(species+"gDNAmap"+str(i+1),
                                 bowtie2_mapping(readfile, "./gDNA/{s}/{s}.{i}.gDNA.sam".format(s=species, i=i+1),
                                                 index_ref_file))

for species in all_species_RNA:
    for tissue in tissues[species]:
        for i, readfile in enumerate(raw_RNA[species][tissue]):
            gwf.target_from_template(species+"_"+tissue+"_"+"RNAmap"+str(i+1),
                                     bowtie2_mapping(readfile, "./RNA/{s}/{s}.{t}.{i}.gDNA.sam".format(s=species, t=tissue, i=i+1),
                                                     index_ref_file))

def samtools_merge(infile1, infile2, outfile):
    inputs = [infile1, infile2]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = """
    source /com/extra/samtools/1.3/load.sh
    samtools merge {} {} {}
    """.format(outfile, infile1, infile2)

    return inputs, outputs, options, spec


for species in all_species_gDNA:
    gwf.target_from_template(species+"gDNAmerge",
                             samtools_merge("./gDNA/{s}/{s}.{i}.gDNA.sam".format(s=species, i=1),
                                            "./gDNA/{s}/{s}.{i}.gDNA.sam".format(s=species, i=2),
                                            "./gDNA/{s}/{s}.gDNA.sam".format(s=species)))

for species in all_species_RNA:
    for tissue in tissues[species]:
        gwf.target_from_template(species+"_"+tissue+"_"+"RNAmerge",
                                 samtools_merge("./RNA/{s}/{s}.{t}.{i}.gDNA.sam".format(s=species, t=tissue, i=1),
                                                "./RNA/{s}/{s}.{t}.{i}.gDNA.sam".format(s=species, t=tissue, i=2),
                                                "./RNA/{s}/{s}.{t}.gDNA.sam".format(s=species, t=tissue)))


def samtools_process(infile, outfile):
    inputs = [infile]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = """
    source /com/extra/samtools/1.3/load.sh

    #samtools view -Sb {infile} -o {infile}.bam
    samtools sort {infile} -o {outfile}
    samtools index {outfile}
    """.format(infile=infile, outfile=outfile)

    return inputs, outputs, options, spec

for species in all_species_gDNA:
    gwf.target_from_template(species+"gDNAsort",
                             samtools_process("./gDNA/{s}/{s}.gDNA.sam".format(s=species),
                                            "./gDNA/{s}.gDNA.s.bam".format(s=species)))

for species in all_species_RNA:
    for tissue in tissues[species]:
        gwf.target_from_template(species+"_"+tissue+"_"+"RNAsort",
                                 samtools_process("./RNA/{s}/{s}.{t}.gDNA.sam".format(s=species, t=tissue),
                                                "./RNA/{s}.{t}.RNA.s.bam".format(s=species, t=tissue)))

def make_protocol_file(in_files, individual, samples, ploidy, filetype, outfile):
    inputs = in_files
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = ''
    for inf, ind, sample, p, f in zip(in_files, individual, samples, ploidy, filetype):
        spec += 'echo "'+"\t".join([ind, p, sample, f, inf])+'" >> '+outfile+"\n"

    print(spec)

    return inputs, outputs, options, spec


Name = {"To": "P1", "Tp": "P2", "TrR": "Ch"}
Ploidy = {"To": "1", "Tp": "1", "TrR": "2"}

bamfiles = []
individual = []
ploidy_levels = []
seq_type = []

#### MAKE SAMPLES LIST
samples = []


translate = {"To": "To", "Tp": "Tp",
             "S9": "TrR", "S10": "TrR"}

samples_map = {"To": [["St1", "R1", "YL1", "OL1"],
                      ["St2", "R2", "YL2", "OL2"],
                      ["St3", "R3", "YL3", "OL3"]],
               "Tp": [["St1", "R1", "YL1", "OL1"],
                      ["St2", "R2", "YL2", "OL2"],
                      ["St3", "R3", "YL3", "OL3"]],
               "S10": [["St1", "R1", "YL1", "OL1"],
                      ["St2", "R2", "YL2", "OL2"],
                      ["St3", "R3", "YL3", "OL3"]],
               "S9": [["F1", "St1", "R1", "L1"],
                      ["F2", "St2", "R2", "L2"],
                      ["F3", "St3", "R3", "L3"]]}

#species_included = {"To": False, "Tp": False, "TrR": False}

for species in [["S9", "S10"], "To", "Tp"]:
    if type(species)==list:
        sample_n = 1
        for specie in species:
            #print(specie, samples_map[specie])
            for current_samples in samples_map[specie]:
                for tissue in current_samples:
                    #print(specie, tissue)
                    bamfiles.append("./RNA/{s}.{t}.RNA.s.bam".format(s=specie, t=tissue))
                    seq_type.append("RNAseq")
                    individual.append(Name[translate[specie]])
                    ploidy_levels.append(Ploidy[translate[specie]])
                    samples.append("sample"+str(sample_n))
                    sample_n += 1
        for specie in [species[0]]:
            bamfiles.append("./gDNA/{s}.gDNA.s.bam".format(s=translate[specie]))
            seq_type.append("gDNA")
            individual.append(Name[translate[specie]])
            ploidy_levels.append(Ploidy[translate[specie]])
            samples.append("sample"+str(sample_n))
            sample_n += 1
    else:
        sample_n = 1
        for current_samples in samples_map[species]:
            for tissue in current_samples:
                bamfiles.append("./RNA/{s}.{t}.RNA.s.bam".format(s=species, t=tissue))
                seq_type.append("RNAseq")
                individual.append(Name[translate[species]])
                ploidy_levels.append(Ploidy[translate[species]])
                samples.append("sample"+str(sample_n))
                sample_n += 1
        bamfiles.append("./gDNA/{s}.gDNA.s.bam".format(s=translate[species]))
        seq_type.append("gDNA")
        individual.append(Name[translate[species]])
        ploidy_levels.append(Ploidy[translate[species]])
        samples.append("sample"+str(sample_n))


#print(bamfiles)
#print(individual)
#print(ploidy_levels)
#print(seq_type)

gwf.target_from_template("HyLiTEprotocolfile",
                         make_protocol_file(bamfiles,
                                            individual,
                                            samples,
                                            ploidy_levels,
                                            seq_type,
                                            "repens_protocol.txt"))

def make_bam_file(in_files,outfile):
    inputs = in_files
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = ""
    for ind in in_files:
        spec += 'echo "'+ind+'" >> '+outfile+"\n"

    return inputs, outputs, options, spec

gwf.target_from_template("mpileup_bamfile",
                         make_bam_file(bamfiles, "bamfiles.txt"))

def mpileup(reference, bamfiles, outfile):
    inputs = [reference, bamfiles]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '24g',
        'account': 'NChain',
        'walltime': '48:00:00'
    }

    spec = '''
    source /com/extra/samtools/1.3/load.sh

    samtools mpileup -BQ 0 -d 1000000 -f {ref} -b {bamfiles} > {outfile}
    '''.format(ref=reference, bamfiles=bamfiles, outfile=outfile)

    return inputs, outputs, options, spec


gwf.target_from_template("mpileup",
                         mpileup(reference_file, "bamfiles.txt",
                                 "repens_combo.pileup"))


def fix_mpileup(pileup_file, outfile):
    inputs = [pileup_file]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = '''
    python mpileupfix.py {pileup} > {out}
    '''.format(pileup=pileup_file, out=outfile)

    return inputs, outputs, options, spec

gwf.target_from_template("mpileupfix",
                         fix_mpileup("repens_combo.pileup", "repens_fixed.pileup"))


def run_HyLiTE(reference, protocol, pileup_file):
    inputs = [reference, protocol, pileup_file]
    outputs = ["HyLiTE_output/HyLiTE_output.expression.txt"]
    options = {
        'cores': 1,
        'memory': '124g',
        'account': 'NChain',
        'walltime': '64:00:00'
    }

    spec = '''
    source activate HyLiTE

    HyLiTE -v -f {proto} -p {pileup} -n HyLiTE_output -r {ref}
    '''.format(ref=reference, proto=protocol, pileup=pileup_file)

    return inputs, outputs, options, spec


def HacknHyLiTE(reference, protocol, pileup_file, segments):
    inputs = [reference, protocol, pileup_file]
    outputs = []
    options = {
        'cores': 1,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    for i in range(segments):
        outputs.append("./slicing_directories/slicing.subset{}.sh".format(i))

    spec = '''
    source activate HyLiTE

    HacknHyLiTE -n {seg} -o slicing_directories --name slicing -f {proto} -p {pileup} --options "-r {ref}"
    '''.format(seg=segments, ref=reference, proto=protocol, pileup=pileup_file)

    return inputs, outputs, options, spec


#gwf.target_from_template("HyLiTE",
#                         run_HyLiTE(reference_file, "repens_protocol.txt", "repens_fixed.pileup"))


n_segments = 40

gwf.target_from_template("HacknHyLiTE",
                         HacknHyLiTE(reference_file, "repens_protocol.txt", "repens_fixed.pileup",
                                     n_segments))

subsets = []
subsets_out = []
for i in range(n_segments):
    subsets.append("./slicing_directories/slicing.subset{}.sh".format(i))
    subsets_out.append("./slicing_directories/subset{}/slicing.subset{}.expression.txt".format(i, i))


def run_subset(subset, output):
    inputs = [subset]
    outputs = [output]
    options = {
        'cores': 1,
        'memory': '24g',
        'account': 'NChain',
        'walltime': '48:00:00'
    }

    spec = '''
    source activate HyLiTE

    bash {subset}
    '''.format(subset=subset)

    return inputs, outputs, options, spec


for i, subset in enumerate(zip(subsets, subsets_out)):
    gwf.target_from_template("Subset"+str(i),
                             run_subset(subset[0], subset[1]))

def mergeHyLiTE(subset_out, directory, reference):
    inputs = subset_out
    outputs = []
    options = {
        'cores': 1,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '48:00:00'
    }

    spec = '''
    source activate HyLiTE

    cd {dir}

    HyLiTE-merge -r {ref}
    '''.format(dir=directory, ref=reference)

    return inputs, outputs, options, spec


gwf.target_from_template("mergeHyLiTE",
                         mergeHyLiTE(subsets_out, "./slicing_directories/", "."+reference_file))
