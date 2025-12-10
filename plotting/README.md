# Plots of simulation and analytic results

## Files

quantities_of_interest.ipynb contains code that plots the quantities relevant to the central questions we are interested in answering/understanding

auxiliary_quantities.ipynb contains code that plots quantities which help us understand why some quantities of interest have inaccurate analytics or behave not as we expected under some cases and parameter regimes

plot_functions.py contains purely analytical functions of the quantities to plot and functions that read and preprocess simulation results

## Quantities of interest

The fixation probability of newly arising mutations

The expected heterozygosity of the population at a given point in time

The relative allelic contribution to the short-term phenotypic change after a shift per unit mutational input (which equals the allelic contribution to the phenotypic variance per unit mutational input)

The proportion of fixations due to adaptation, $\alpha$, in the McDonald-Kreitman test, with and without allele frequency cutoff

## Auxiliary quantities

The average expected change in allele frequency

The variance of change in allele frequency

The total number of shifts that a fixed mutation experiences

The excess number of aligned shifts than opposing shifts that a fixed mutation experiences
