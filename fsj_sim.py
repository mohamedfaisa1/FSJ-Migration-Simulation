import numpy as np
import msprime as mp
import tskit as tsk
import argparse
import os


# function that returns the demography based on different migration schemes and a set migration rate
def setupDemography(flag: int = 0, m=0.05):
    """
    flag = 0 -> no migration
    flag = 1 -> unidirectional
    flag = 2 -> reciprocal
    flag = 3 -> directional
    m -> migration rate (generally between LARGEPOP & SMALLPOP, will be empirically scaled for LARGEPOP & MEDIUMPOP)
    """
    demography = mp.Demography()

    # Basic Setup (Present Day)
    demography.add_population(
        name="ANC",
        description="Ancestral Population",
        initial_size=15000,
    )
    demography.add_population(
        name="LARGEPOP",
        initial_size=1000
    )
    demography.add_population(
        name="MEDIUMPOP",
        initial_size=300
    )
    demography.add_population(
        name="SMALLPOP",
        initial_size=50
    )

    # 1912 Railroad 15 Generations Ago
    demography.add_population_parameters_change(time=15, population="LARGEPOP", initial_size=2500)
    demography.add_population_parameters_change(time=15, population="Medium", initial_size=1000)
    demography.add_population_parameters_change(time=15, population="SMALLPOP", initial_size=1500)

    # European Colonization End 40 Generations Ago
    demography.add_population_parameters_change(time=40, population="LARGEPOP", initial_size=4000)
    demography.add_population_parameters_change(time=40, population="MEDIUMPOP", initial_size=2000)
    demography.add_population_parameters_change(time=40, population="SMALLPOP", initial_size=3000)

    # European Colonization Start 80 Generations Ago
    demography.add_population_parameters_change(time=80, population="LARGEPOP", initial_size=10000)
    demography.add_population_parameters_change(time=80, population="MEDIUMPOP", initial_size=6000)
    demography.add_population_parameters_change(time=80, population="SMALLPOP", initial_size=7000)

    # Ancestral Fragmentation 2000 Generations Ago
    demography.add_population_split(time=2000, derived=["LARGEPOP", "MEDIUMPOP", "SMALLPOP"], ancestral="ANC")
    demography.add_population_parameters_change(time=2000, population="ANC", initial_size=15000)

    # Ancestral Population Decline 8000 Generations Ago
    demography.add_population_parameters_change(time=8000, population="ANC", initial_size=45000)

    ##  Migration

    # case (default): no migration
    migration_matrix = np.zeros(16).reshape(4, 4)

    # case: unidirectional
    if flag == 1:
        migration_matrix[3, 1] = m
        migration_matrix[2, 1] = m * 0.5
    # case: reciprocal
    elif flag == 2:
        migration_matrix[3, 1] = m
        migration_matrix[1, 3] = m
        migration_matrix[2, 1] = m * 0.5
        migration_matrix[1, 2] = m * 0.5
    # case: directional
    elif flag == 3:
        migration_matrix[3, 1] = m
        migration_matrix[2, 1] = m * 0.5
        migration_matrix[1, 3] = m * 0.5
        migration_matrix[1, 2] = m * 0.25

    demography.migration_matrix = migration_matrix

    demography.sort_events()  # ensures all events are sorted (to avoid bugs)

    return demography


# function to calculate Watterson's Theta from segregating sites
def wattersonTheta(sim, L, n):
    S = tsk.TreeSequence.segregating_sites(sim, mode="site")
    hs = np.sum(1/np.arange(1, n))
    theta = S / (hs * L)
    return theta, S


# just a function to output population statistics for each population
def outputPopulationStatistics(sim, N: int, n: int, L: int, s:int, name:str, seed=1989):
    """
    :param sims: this is a generator object containing N replicates of simulations
    :return: returns 3 vectors, each of length N, which are Tajima's D, PI (genetic) diversity,
    and Watterson's Theta
    """

    tajD = np.zeros(3)
    piDiv = np.zeros(3)
    watterson = np.zeros(3)
    seg_sites = np.zeros(3)
    sfs_data = np.zeros((3, (2*n)))
    hs = np.sum(1/np.arange(1,n))

    mm = mp.sim_mutations(sim, 3e-9, random_seed=seed, model=mp.JC69())

    LARGEPOP = mm.samples(population=1)
    MEDIUMPOP = mm.samples(population=2)
    SMALLPOP = mm.samples(population=3)

    LARGEPOP_sims = mm.subset(LARGEPOP)
    MEDIUMPOP_sims = mm.subset(MEDIUMPOP)
    SMALLPOP_sims = mm.subset(SMALLPOP)

    sfs_data[0] = LARGEPOP_sims.allele_frequency_spectrum(mode="site", polarised=True,
                                                    span_normalise=False)[1:]
    sfs_data[1] = MEDIUMPOP_sims.allele_frequency_spectrum(mode="site", polarised=True,
                                                    span_normalise=False)[1:]
    sfs_data[2] = SMALLPOP_sims.allele_frequency_spectrum(mode="site", polarised=True,
                                                    span_normalise=False)[1:]

    piDiv[0] = LARGEPOP_sims.diversity(mode="site")
    piDiv[1] = MEDIUMPOP_sims.diversity(mode="site")
    piDiv[2] = SMALLPOP_sims.diversity(mode="site")

    seg_sites[0] = LARGEPOP_sims.segregating_sites(mode="site")
    seg_sites[1] = MEDIUMPOP_sims.segregating_sites(mode="site")
    seg_sites[2] = SMALLPOP_sims.segregating_sites(mode="site")

    watterson[0] = seg_sites[0] / (hs*L)
    watterson[1] = seg_sites[1] / (hs*L)
    watterson[2] = seg_sites[2] / (hs*L)

    tajD[0] = LARGEPOP_sims.Tajimas_D(mode="site")
    tajD[1] = MEDIUMPOP_sims.Tajimas_D(mode="site")
    tajD[2] = SMALLPOP_sims.Tajimas_D(mode="site")

    Fst_OA = mm.Fst([LARGEPOP, MEDIUMPOP])
    Fst_OC = mm.Fst([LARGEPOP, SMALLPOP])
    Fst_AC = mm.Fst([MEDIUMPOP, SMALLPOP])

    LARGEPOP_f = open("LARGEPOP__"+ name + str(s + 1) + ".vcf", 'w')
    LARGEPOP_sims.write_vcf(output=LARGEPOP_f)
    LARGEPOP_f.close()

    MEDIUMPOP_f = open("MEDIUMPOP__" + name + str(s + 1) + ".vcf", 'w')
    MEDIUMPOP_sims.write_vcf(output=MEDIUMPOP_f)
    MEDIUMPOP_f.close()

    SMALLPOP_f = open("SMALLPOP__" + name + str(s + 1) + ".vcf", 'w')
    SMALLPOP_sims.write_vcf(output=SMALLPOP_f)
    SMALLPOP_f.close()

    return tajD, watterson, piDiv, seg_sites, sfs_data, Fst_OA, Fst_OC, Fst_AC



# function to set up sampling scheme
"""
sample_history = False -> sample from present-day (time = 0)
sample_history = True -> sample throughout the past in a predetermined set of intervals
sample_total = False -> sample subset of individuals from samples at present-day (time = 0)
sample_total = True -> sample all individuals form samples
"""


def setupSamplingSchemes(sample_history: bool = False, sample_total: bool = False):
    n = np.zeros(3, dtype=np.int32)  # list of all sample sizes
    sample_schemes = []

    n[0] = 50
    n[1] = 50
    n[2] = 50
    sample_LARGEPOP = [mp.SampleSet(50, "LARGEPOP", 0)]
    sample_MEDIUMPOP = [mp.SampleSet(50, "MEDIUMPOP", 0)]
    sample_SMALLPOP = [mp.SampleSet(50, "SMALLPOP", 0)]
    sample_schemes.extend(sample_LARGEPOP)
    sample_schemes.extend(sample_MEDIUMPOP)
    sample_schemes.extend(sample_SMALLPOP)
    return sample_schemes, n



def main():
    # For Command Line Arguments
    parser = argparse.ArgumentParser(
        description="Simulate Florida scrub-jay demographic history with varying migration schemes" +
                    " and calculate population genetic summary statistics.")
    parser.add_argument("-replicates", action="store", dest="replicates", type=int,
                        default=100, help="Number of replicate simulations to perform.")
    parser.add_argument("-recombination_rate", action="store", dest="recombination_rate", type=float,
                        default=1.5e-8, help="Recombination rate.")
    parser.add_argument("-migration", action="store", dest="m", type=float,
                        default=0.05, help="Migration rate.")
    parser.add_argument("-mutations", action="store", dest="mut", type=bool,
                        default=True, help="Simulate Jukes-Cantor 1969 Model of Mutations.")
    parser.add_argument("-L", action="store", dest="L", type=float,
                        default=1e6, help="Length of genomic segment to simulate.")
    parser.add_argument("-end", action="store", dest="end", type=int,
                        default=None, help="End time (generations) of simulation.")
    parser.add_argument("-seed", action="store", dest="seed", type=int,
                        default=19, help="Random seed.")
    parser.add_argument("-sample_history", action="store", dest="sample_history", type=bool,
                        default=False,
                        help="If True, the sampling scheme will sample throughout history of the simulation." +
                             " [Note: If True, sample_total should be set to False.]")
    parser.add_argument("-sample_total", action="store", dest="sample_total", type=bool,
                        default=True,
                        help="If True, the sampling scheme will sample the whole population " +
                             "from the present-day (generation 0). If False, it will sample down.")
    args = parser.parse_args()

    # Useful Constants
    modes = ["site", "branch"]  # useful for indicating which mode to calculate population genetic statistics
    replicates = args.replicates  # number of replicates for each simulation
    rec = args.recombination_rate  # genomic segment recombination rate
    m = args.m # migration rate
    L = args.L  # length of genomic segment
    end = args.end  # time to end simulation
    seed = args.seed  # random seed
    sample_history = args.sample_history # whether to sample historically
    sample_total = args.sample_total # whether to sample total or down
    mutations = args.mut # whether to simulate mutations
    mu = 3e-9 # mutation rate
    JC69 = mp.JC69() # Jukes-Cantor Model 1969
    populations = ["LARGEPOP", "MEDIUMPOP", "SMALLPOP"]

    sample_schemes, n = setupSamplingSchemes()

    # Initialize Demography Scenarios
    demography_0 = setupDemography(flag=0, m=0) # no migration scheme
    demography_A = setupDemography(flag=1, m=m) # unidirectional migration scheme
    demography_B = setupDemography(flag=2, m=m) # reciprocal migration scheme
    demography_C = setupDemography(flag=3, m=m) # directional migration scheme

    # Sample Different from Demographic Schemes
    sim_NM = mp.sim_ancestry(samples=sample_schemes, demography=demography_0,
                             sequence_length=L, recombination_rate=rec, end_time=end,
                             num_replicates=replicates, random_seed=seed)
    sim_UD = mp.sim_ancestry(samples=sample_schemes, demography=demography_A,
                             sequence_length=L, recombination_rate=rec, end_time=end,
                             num_replicates=replicates, random_seed=seed)
    sim_RC = mp.sim_ancestry(samples=sample_schemes, demography=demography_B,
                             sequence_length=L, recombination_rate=rec, end_time=end,
                             num_replicates=replicates, random_seed=seed)
    sim_DM = mp.sim_ancestry(samples=sample_schemes, demography=demography_C,
                             sequence_length=L, recombination_rate=rec, end_time=end,
                             num_replicates=replicates, random_seed=seed)

    # FINAL RESULTS via Migration Scheme

    ## Setting up Population Statistic Arrays
    NM_Fst_OA = np.zeros(replicates)
    NM_Fst_OC = np.zeros(replicates)
    NM_Fst_AC = np.zeros(replicates)

    UD_Fst_OA = np.zeros(replicates)
    UD_Fst_OC = np.zeros(replicates)
    UD_Fst_AC = np.zeros(replicates)

    RC_Fst_OA = np.zeros(replicates)
    RC_Fst_OC = np.zeros(replicates)
    RC_Fst_AC = np.zeros(replicates)

    DM_Fst_OA = np.zeros(replicates)
    DM_Fst_OC = np.zeros(replicates)
    DM_Fst_AC = np.zeros(replicates)

    sfs_0 = np.zeros((3, 2 * n[0]))
    sfs_A = np.zeros((3, 2 * n[1]))
    sfs_B = np.zeros((3, 2 * n[2]))
    sfs_C = np.zeros((3, 2 * n[2]))

    tajD_0 = np.zeros(shape=(3, replicates))
    pi_0 = np.zeros(shape=(3, replicates))
    watterson_0 = np.zeros(shape=(3, replicates))
    seg_0 = np.zeros(shape=(3, replicates))

    tajD_A = np.zeros(shape=(3, replicates))
    pi_A = np.zeros(shape=(3, replicates))
    watterson_A = np.zeros(shape=(3, replicates))
    seg_A = np.zeros(shape=(3, replicates))

    tajD_B = np.zeros(shape=(3, replicates))
    pi_B = np.zeros(shape=(3, replicates))
    watterson_B = np.zeros(shape=(3, replicates))
    seg_B = np.zeros(shape=(3, replicates))

    tajD_C = np.zeros(shape=(3, replicates))
    pi_C = np.zeros(shape=(3, replicates))
    watterson_C = np.zeros(shape=(3, replicates))
    seg_C = np.zeros(shape=(3, replicates))

    # get results

    for s, sim in enumerate(sim_NM):
        (tajD_0[:,s], watterson_0[:,s], pi_0[:,s], seg_0[:,s],
           sfs_h, NM_Fst_OA[s], NM_Fst_OC[s], NM_Fst_AC[s]) = outputPopulationStatistics(sim, replicates, n[0],
                                                              L, s, "NM__muts__" + str(m) + "_")
        sfs_0 += sfs_h

    for s, sim in enumerate(sim_UD):
        (tajD_A[:,s], watterson_A[:,s], pi_A[:,s], seg_A[:,s],
         sfs_h, UD_Fst_OA[s], UD_Fst_OC[s], UD_Fst_AC[s]) = outputPopulationStatistics(sim, replicates, n[0],
                                                                L, s, "UD__muts__" + str(m) + "_")
        sfs_A += sfs_h

    for s, sim in enumerate(sim_RC):
        (tajD_B[:,s], watterson_B[:,s], pi_B[:,s], seg_B[:,s],
         sfs_h, RC_Fst_OA[s], RC_Fst_OC[s], RC_Fst_AC[s]) = outputPopulationStatistics(sim, replicates, n[0],
                                                                L, s, "RC__muts__" + str(m) + "_")
        sfs_B += sfs_h

    for s, sim in enumerate(sim_DM):
        (tajD_C[:,s], watterson_C[:,s], pi_C[:,s], seg_C[:,s],
         sfs_h, DM_Fst_OA[s], DM_Fst_OC[s], DM_Fst_AC[s]) = outputPopulationStatistics(sim, replicates, n[0],
                                                                L, s, "DM__muts__" + str(m) + "_")
        sfs_C += sfs_h


    # Tajima's D
    np.save("tajimaD_no_mig_" + str(m), tajD_0)
    np.save("tajimaD_unidirectional_" + str(m), tajD_A)
    np.save("tajimaD_reciprocal_" + str(m), tajD_B)
    np.save("tajimaD_directional_" + str(m), tajD_C)

    # Pi (Genetic Diversity)
    np.save("pi_no_mig_" + str(m), pi_0)
    np.save("pi_unidirectional_" + str(m), pi_A)
    np.save("pi_reciprocal_" + str(m), pi_B)
    np.save("pi_directional_" + str(m), pi_C)

    # Watterson's Theta
    np.save("wt_no_mig_" + str(m), watterson_0)
    np.save("wt_unidirectional_" + str(m), watterson_A)
    np.save("wt_reciprocal_" + str(m), watterson_B)
    np.save("wt_directional_" + str(m), watterson_C)

    # Seg Sites
    np.save("seg_no_mig_" + str(m), seg_0)
    np.save("seg_unidirectional_" + str(m), seg_A)
    np.save("seg_reciprocal_" + str(m), seg_B)
    np.save("seg_directional_" + str(m), seg_C)

    # SFS
    np.save("sfs_nm_" + str(m), sfs_0)
    np.save("sfs_ud_" + str(m), sfs_A)
    np.save("sfs_rc_" + str(m), sfs_B)
    np.save("sfs_dm_" + str(m), sfs_C)

    # Fst
    np.save("OA_no_mig_" + str(m), NM_Fst_OA)
    np.save("OC_no_mig_" + str(m), NM_Fst_OC)
    np.save("AC_no_mig_" + str(m), NM_Fst_AC)

    np.save("OA_unid_" + str(m), UD_Fst_OA)
    np.save("OC_unid_" + str(m), UD_Fst_OC)
    np.save("AC_unid_" + str(m), UD_Fst_AC)

    np.save("OA_recip_" + str(m), RC_Fst_OA)
    np.save("OC_recip_" + str(m), RC_Fst_OC)
    np.save("AC_recip_" + str(m), RC_Fst_AC)

    np.save("OA_dirmig_" + str(m), RC_Fst_OA)
    np.save("OC_dirmig_" + str(m), RC_Fst_OC)
    np.save("AC_dirmig_" + str(m), RC_Fst_AC)



if __name__ == "__main__":
    main()
