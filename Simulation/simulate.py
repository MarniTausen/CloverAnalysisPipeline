import msprime
from optparse import OptionParser

def help():
    print "====== Simulation through msprime with a simple growing population ====="
    print ""
    print "-o <filename.vcf>                           Output vcf file"
    print "-N int                                      Effective population size"
    print "-m float                                    Mutation rate"
    print "-r float                                    Recombination rate"
    print "-l int                                      Genome length"
    print "-g float                                    Growth rate"
    print "-s float                                    Sample size"
    print "-t int                                      Time"
    print "--print                                     Printing debug info"
    exit()

if __name__=="__main__":

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename.vcf>")
    parser.add_option('-N', type="float", nargs=1, dest="Popsize", help="Effective population size")
    parser.add_option('-m', type="float", nargs=1, dest="mutrate", help="Mutation rate")
    parser.add_option('-r', type="float", nargs=1, dest="recomrate", help="Recombination rate")
    parser.add_option('-l', type="float", nargs=1, dest="length", help="Length")
    parser.add_option('-g', type="float", nargs=1, dest="growth", help="growth rate")
    parser.add_option('-s', type="int", nargs=1, dest="samples", help="Sample size")
    parser.add_option('-t', type="float", nargs=1, dest="time", help="Time")
    parser.add_option('-d', type="int", nargs=1, dest="timeoff", help="Time of population decreasion")
    parser.add_option('-b', type="int", nargs=1, dest="bottleneck", help="Bottleneck size")
    parser.add_option('--print', action="store_true", dest="prints", help="Print debug")
    parser.add_option('--start', action="store_true", dest="start", help="Trigger Population start with bottleneck")
    parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
    options, args = parser.parse_args()

    if options.help!=None:
        help()

    if options.Popsize!=None:
        Ne = options.Popsize
    else:
        Ne = 10000
    if options.mutrate!=None:
        mu = options.mutrate
    else:
        print "No mutation rate given. This will produce no variants."
        print "use -m to add a mutation rate"
        mu = None
    if options.recomrate!=None:
        r = options.recomrate
    else:
        r = None
    if options.length!=None:
        L = options.length
    else:
        L = None
    if options.growth!=None:
        g = options.growth
    else:
        g = 0
    if options.samples!=None:
        ss = options.samples
    else:
        ss = 1
    if options.time!=None:
        t = options.time
    else:
        t = 1
    if options.timeoff!=None:
        d = options.timeoff
    else:
        d = None
    if options.bottleneck!=None:
        bottleneck_size = options.bottleneck
    else:
        bottleneck_size = None


    population_configurations = [msprime.PopulationConfiguration(sample_size=ss,
                                                                 initial_size=Ne,
                                                                 growth_rate=0)]

    if d==None and bottleneck_size==None:
        demographic_events = [msprime.PopulationParametersChange(t, growth_rate=0)]
    elif d!=None and bottleneck_size==None:
        demographic_events = [msprime.PopulationParametersChange(t-d, growth_rate=g),
                              msprime.PopulationParametersChange(t, growth_rate=0)]
    elif d==None and bottleneck_size!=None:
        if options.start!=None:
            demographic_events = [msprime.PopulationParametersChange(t, initial_size=bottleneck_size)]
        else:
            demographic_events = [msprime.PopulationParametersChange(t, initial_size=bottleneck_size),
                                  msprime.PopulationParametersChange(t+1, initial_size=Ne)]
    else:
        raise """
Warning, using -d and -b together are not implemented. For instanstaneous bottlenecks use -b.
Otherwise use -d and -g to get the bottleneck size you want."""

    if options.prints!=None:
        dd = msprime.DemographyDebugger(Ne=Ne,
                                        population_configurations=population_configurations,
                                        demographic_events=demographic_events)
        dd.print_history()

    TS = msprime.simulate(Ne=Ne, length=L, recombination_rate=r, mutation_rate=mu,
                          population_configurations=population_configurations,
                          demographic_events=demographic_events)

    if options.output!=None:
        out = options.output
    else:
        out = "output.vcf"

    with open(out, "w") as vcf_file:
        TS.write_vcf(vcf_file, 2)


