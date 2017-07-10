from gwf import Workflow

gwf = Workflow()

def bwa_map(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '4g',
        'cores': '8',
        'walltime': '240:00:00',
        'account': 'NChain'
    }

    spec = "./map.sh "+" ".join(infiles+outfiles)

    return options, spec

gwf.target("NCL08") << bwa_map(["NCL-08/FCHGMFMBBXX-wHAXPI040514-50_L8_1.fq.gz",
                                "NCL-08/FCHGMFMBBXX-wHAXPI040514-50_L8_2.fq.gz",
                                "NCL-08/FCHGMFNBBXX-wHAXPI040514-50_L2_1.fq.gz",
                                "NCL-08/FCHGMFNBBXX-wHAXPI040514-50_L2_2.fq.gz"],
                               ["ncl-08L8.bam", "ncl-08L2.bam"])

gwf.target("NCL09") << bwa_map(["NCL-09/FCHGMFMBBXX-wHAXPI040513-52_L8_1.fq.gz",
                                "NCL-09/FCHGMFMBBXX-wHAXPI040513-52_L8_2.fq.gz",
                                "NCL-09/FCHGMFNBBXX-wHAXPI040513-52_L2_1.fq.gz",
                                "NCL-09/FCHGMFNBBXX-wHAXPI040513-52_L2_2.fq.gz"],
                               ["ncl-09L8.bam", "ncl-09L2.bam"])

gwf.target("NCL10") << bwa_map(["NCL-10/FCHGMFMBBXX-wHAXPI040512-53_L8_1.fq.gz",
                                "NCL-10/FCHGMFMBBXX-wHAXPI040512-53_L8_2.fq.gz",
                                "NCL-10/FCHGMFNBBXX-wHAXPI040512-53_L2_1.fq.gz",
                                "NCL-10/FCHGMFNBBXX-wHAXPI040512-53_L2_2.fq.gz"],
                               ["ncl-10L8.bam", "ncl-10L2.bam"])

gwf.target("NCL12") << bwa_map(["NCL-12/FCHGMFMBBXX-wHAXPI040511-54_L8_1.fq.gz",
                                "NCL-12/FCHGMFMBBXX-wHAXPI040511-54_L8_2.fq.gz",
                                "NCL-12/FCHGMFNBBXX-wHAXPI040511-54_L2_1.fq.gz",
                                "NCL-12/FCHGMFNBBXX-wHAXPI040511-54_L2_2.fq.gz"],
                               ["ncl-12L8.bam", "ncl-12L2.bam"])


def merge_reads(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '4g',
        'walltime': '144:00:00',
        'account': 'NChain'
    }

    spec = "./merge.sh "+" ".join(infiles+outfiles)

    return options, spec

gwf.target("Merge_NCL08") << merge_reads(["ncl-08L8.bam", "ncl-08L2.bam"],
                                         ["ncl-08.u.bam"])

gwf.target("Merge_NCL09") << merge_reads(["ncl-09L8.bam", "ncl-09L2.bam"],
                                         ["ncl-09.u.bam"])

gwf.target("Merge_NCL10") << merge_reads(["ncl-10L8.bam", "ncl-10L2.bam"],
                                         ["ncl-10.u.bam"])

gwf.target("Merge_NCL12") << merge_reads(["ncl-12L8.bam", "ncl-12L2.bam"],
                                         ["ncl-12.u.bam"])

def sort_reads(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '4g',
        'walltime': '144:00:00',
        'account': 'NChain'
    }

    spec = "./sort.sh "+" ".join(infiles+outfiles)

    return options, spec

gwf.target("Sort_NCL08") << sort_reads(["ncl-08.u.bam"],
                                       ["ncl-08.bam"])

gwf.target("Sort_NCL09") << sort_reads(["ncl-09.u.bam"],
                                       ["ncl-09.bam"])

gwf.target("Sort_NCL10") << sort_reads(["ncl-10.u.bam"],
                                       ["ncl-10.bam"])

gwf.target("Sort_NCL12") << sort_reads(["ncl-12.u.bam"],
                                       ["ncl-12.bam"])


def indexbam(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '1g',
        'walltime': '144:00:00',
        'account': 'NChain'
    }

    spec = "./indexbam.sh "+" ".join(infiles+outfiles)

    return options, spec

gwf.target("Index_NCL08") << indexbam(["ncl-08.bam"],
                                      ["ncl-08.bai"])

gwf.target("Index_NCL09") << indexbam(["ncl-09.bam"],
                                      ["ncl-09.bai"])

gwf.target("Index_NCL10") << indexbam(["ncl-10.bam"],
                                      ["ncl-10.bai"])

gwf.target("Index_NCL12") << indexbam(["ncl-12.bam"],
                                      ["ncl-12.bai"])

def bamToDiploid(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '1g',
        'walltime': '72:00:00',
        'account': 'NChain'
    }

    spec = "./bam2diploid.sh "+" ".join(infiles+outfiles)

    return options, spec

gwf.target("diploid_NCL08") << bamToDiploid(["ncl-08.bam"],
                                            ["ncl-08.fq.gz"])

gwf.target("diploid_NCL09") << bamToDiploid(["ncl-09.bam"],
                                            ["ncl-09.fq.gz"])

gwf.target("diploid_NCL10") << bamToDiploid(["ncl-10.bam"],
                                            ["ncl-10.fq.gz"])

gwf.target("diploid_NCL12") << bamToDiploid(["ncl-12.bam"],
                                            ["ncl-12.fq.gz"])

def PSMC(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '1g',
        'walltime': '24:00:00',
        'account': 'NChain'
    }

    spec = "./psmc.sh "+" ".join(infiles)
    
    return options, spec

gwf.target("psmc_NCL08") << PSMC(["ncl-08.fq.gz"],
                                 ["ncl-08.fq.psmc", "ncl-08.fq.psmcfa"])

gwf.target("psmc_NCL09") << PSMC(["ncl-09.fq.gz"],
                                 ["ncl-09.fq.psmc", "ncl-09.fq.psmcfa"])

gwf.target("psmc_NCL10") << PSMC(["ncl-10.fq.gz"],
                                 ["ncl-10.fq.psmc", "ncl-10.fq.psmcfa"])

gwf.target("psmc_NCL12") << PSMC(["ncl-12.fq.gz"],
                                 ["ncl-12.fq.psmc", "ncl-12.fq.psmcfa"])


def CombinePSMC(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '1g',
        'walltime': '1:00:00',
        'account': 'NChain'
    }

    spec = "cat "+" ".join(infiles)+">> "+outfiles[0]
    
    return options, spec


gwf.target("combined") << CombinePSMC(["ncl-08.fq.psmc",
                                       "ncl-09.fq.psmc",
                                       "ncl-10.fq.psmc",
                                       "ncl-12.fq.psmc"],
                                     ["combined.psmc"])

def PlotPSMC(infiles, outfiles):
    options = {
        'inputs': infiles,
        'outputs': outfiles,
        'memory': '1g',
        'walltime': '24:00:00',
        'account': 'NChain'
    }

    spec = "./plotpsmc.sh "+" ".join(infiles)

    return options, spec

gwf.target("plotPSMC") << PlotPSMC(["combined.psmc"],
                                   ["combined.eps"])

