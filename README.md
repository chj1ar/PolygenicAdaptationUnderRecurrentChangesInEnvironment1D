# Polygenic Adaptation under Recurrent Changes in Environment
### This is the code for the simulations from ???


# Overview

The purpose of the code is to simulate a constant size diploid population that is at steady-state under stabilizing selection and recurrent shifts in fitness optimum. For further details about the scenario and the theory behind the simulations please see ???

The current folder contains code simulating the populations evolving, code running the simulations, and code processing the simulation results. The simulation code records the state of a simulation from time to time, in case the simulation crashes, and resumes the simulation from the recorded state if so.

These programs can be run using command line.

## Simulations

populations\_functions\_record\_and\_resume.py contains the classes to simulate the population.

mutation\_functions.py has container classes for mutations arising in the population.

## To run simulations

Use populations\_argparse\_recurrent\_shifts\_record.py and populations\_argparse\_recurrent\_shifts\_no\_efs\_bins.py to run simulations with and without effect size bins, respectively. The programs can be run on the command line (see population_recurrent_shifts.sh as an example) and take the following parameters:

-N population size

-Vs squared width of fitness function. If you choose Vs to be negative it will automatically be set to 2N (default=-1)

-U mutation rate per haplotype genome per generation

-D\_s0, --shift\_s0 the shift in units of sigma\_0 (note that (sigma_0)<sup>2</sup> = V<sub>A</sub> in the manuscript), which is the expected standard deviation of the phenotype distribution at steady-state under stabilizing selection

-E2Ns expected value of the gamma distribution of scaled selection coefficients (which are the squared phenotypic effects in units of Vs/(2N))

-V2Ns variance of the gamma distribution of scaled selection coefficients. If you choose Vs to be negative it will automatically be set to E2Ns<sup>2</sup>, making the distribution into an exponential distribution with expected value E2Ns (default=-1)

-nE number of effect size bins over which we average the fixation probability

-nF number of frequency bins (currently unused)

-bTN, --burn_time_N The burntime in units of population size (recommend 10)

--rate_of_shift the rate of shifts, i.e., the probability that a shift occurs at each given generation

-nnM, --number_new_mutations number of new mutations to calculate the fixation probability

--complete-jumps specifying this means that there is no shift when the mean phenotype approaches the new optimum after a shift (short-term phenotypic adaptation), to test whether the discrepancies between all allele simulation results and diffusion approximation are because of shifts during the short-term phenotypic adaptation

-sigma_0_del specify this if you want to run all allele simulations using the empirical V<sub>A</sub> from some previous iteration of this simulation

-sD, --save_directory The directory in which you want your simulation results to be saved

--directories_for_empirical_V_A The simulation result directories where the empirical V_A is used in this simulation. We input only the directory with the last array index of the directories, but will use all the directories

## Simulation results

The file 'identifiers.txt' inside the save_directory records the choices of simulation parameters.

The file 'fixation_probability' (with effect size bins) and 'fixations_and_extinctions' (without effect size bins) inside the save_directory records statistics in order to calculate the fixation probabilities of newly arising mutations, the expected heterozygosity of the population at a given point in time, and the relative allelic contribution to the short-term phenotypic change after a shift per unit mutational input, and the proportion of fixations due to adaptation, $\alpha$, in the McDonald-Kreitman test.

## To process simulation results

plotting/plot_functions.py contains functions that calculate all the relevant quantities from simulation results. For quantities whose calculations take long, we have the following scripts that run the functions and save the results (of the quantities). The specific quantities and parameters are elaborated in the scripts.

analytic_fixation_probability_shift_size_gamma.py

analytic_heterozygosity.py

analytic_heterozygosity_truncating_rare_alleles.py

contribution_to_change.py

contribution_to_change_across_rates_of_shifts.py

contribution_to_change_across_shift_sizes.py

average_E_Delta_x.py

average_E_Delta_x_across_rates_of_shifts.py

average_E_Delta_x_across_shift_sizes.py

number_of_shifts_a_fixed_mutation_experiences.py

excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences.py

analytic_V_A.py

## To resume a simulation

If a simulation crashes, then we can resume it with populations_argparse_recurrent_shifts_resume.py, which first loads the recorded state and then runs the simulation from the recorded state. The program can be run on the command line and takes the following parameters:

-iF, --intermediate_Folder The directory in which the recorded state is saved. This directory will also save your simulation results.
