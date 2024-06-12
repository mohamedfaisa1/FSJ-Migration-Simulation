# FSJ-Migration-Simulation
Simulating the Evolution of Neutral Genetic Diversity in the Federally Threatened Florida Scrub-Jay via Varying Migratory Schemes and Rates

## Usage

`fsj_sim.py` outputs NumPy files (32 in total) for each of the summary statistics calculated.

```
python3 fsj_sim.py [-h] [-replicates N] [-recombination_rate rho] 
                        [-migration m] [-mutations mu_bool] [-L L]
                        [-end END] [-seed seed] [-sample_history sh_bool]
                        [-sample_total st_bool]
```

Command-line arguments:
- `-replicates`: the number of replicates of simulations to perform | default = 100 
- `-recombination_rate`: the recombination rate | default: $\rho$ = 1.5e-8
- `-migration`: the migration rate | default: m = 0.05
- `-mutations`: a boolean varibale; indicates whether or not to simualte Jukes-Cantor 1969 style mutations with $\mu$ = 3e-9 | default = True
- `-L`: the length of the macrochromosome to simualte | default = 1e7
- `-end`: the generation for which the simulation should end | default = None
- `-seed`: the random seed for simualtions | default = 1
- `-sample_history`: a boolean variable; if True, the sampling scheme will sample throughout history of the simulation.
                    [Note: If True, sample_total should be set to False.] | default = False
- `-sample_total`: a boolean variable; if True, the sampling scheme will sample the whole population from the present-day (generation 0).
                    if False, it will sample down | default = True


`vcf_calcs.py` computes the number of heterozygous and homozygous sites per population per migratory scheme per individual. Results in a total of 24 NumPy files.
```
python3 vcf_calcs.py [-h] [-dir DIR] [-n n] 
```

Command-line arguments:
- `-dir`: the directory of the vcf files (separated by migratory scheme and population)
- `-n`: the number of samples in the vcf file | default = 50
