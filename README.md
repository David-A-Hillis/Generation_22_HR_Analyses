# High-Runner selection experiment allelic comparisons between generations 22 and 61

This code is based those used in:
	David A Hillis, Liran Yadgary, George M. Weinstock, Fernando Pardo-Manuel de Villena, Daniel Pomp, 
	and Theodore Garland Jr. “Large Changes in Detected Selection Signatures After a Selection Limit 
	in Mice Bred for Voluntary Wheel-Running Behavior.” PLOS ONE, 2024.

Data can be found either in the repositories described by the publications or at: 
	https://drive.google.com/drive/u/1/folders/17xSxKyeow9FeMtOJFRn33ru976EUIgjD
These data are allele frequency data for "WRF.R" under the "Statistical_and_Power_Analyses".
Data included in these files are chromosome (column 1), base pair position (column 2), allele frequencies 
for control lines (columns 3 through 6), and allele frequencies for High-Runner lines (columns 7 through 
10).

# Simulations
Included here is code for conducting simulations of changes in allele frequency in unlinked loci subject
in response to within family selection. These simulations consider family size, total pairings, both 
genetic and environmental variance, different effect sizes of loci, and more. These simulations are 
constructed to make the running outputs, heritability, and seasonal variation as realistic as possible.

Epistasis with loci affecting the constraint or season can be changed, but as default do not. 

As default, 5 constrained and 5 unconstrained simulations are performed (total runtime per simulation can 
be 6 hours due to tracking 2096 loci).

Included are additional scripts for calculating power, selection differential, and heritability from the 
results. Graphical representations are also included in the scripts. Modifications of some of the 
simulation parameters will create the need for modifications in some of these additional scripts.

## Simulations Future Plans

1 - Include additional analytical code from this manuscript

2 - Update simulation code so that it is modular and not a bunch of loops

3 - Update simulations to generate fewer intermediary files

4 - Update analytical code for the simulations to better incorporate different simulation options 

# Statistical and Power Analyses
The code in this folder will conduct the windowed, regularized F-test described in Hillis et al. 2024.

## Statistical and Power Analyses - Future Plans

1 - Add modifications for parallelization

2 - Add code for power analyses