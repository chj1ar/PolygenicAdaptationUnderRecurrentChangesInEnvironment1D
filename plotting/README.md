# Plots of simulation and analytic results

## Files

plot_functions.py contains purely analytical functions of the quantities to plot and functions that read and preprocess simulation results

each other script is for one distinct type of figure. One can find the documentation at the start of each script.

The quantities that the scripts plot are listed below.

## Quantities of interest

The fixation probability of newly arising mutations

The expected heterozygosity of the population at a given point in time

The relative allelic contribution to the short-term phenotypic change after a shift per unit mutational input (which equals the allelic contribution to the phenotypic variance per unit mutational input)

The proportion of fixations due to adaptation, $\alpha$, in the McDonald-Kreitman test, with and without allele frequency cutoff

## Phenotypic dynamics

The distance of the (population) phenotypic mean to the fitness optimum, $D$

The phenotypic variance, $V_A$

The third central moment of the phenotype, directly measuring the phenotypic skewness, $\mu_3$

## Auxiliary quantities

The average expected change in allele frequency

The variance of change in allele frequency

The total number of shifts that a fixed mutation experiences

The excess number of aligned shifts than opposing shifts that a fixed mutation experiences
