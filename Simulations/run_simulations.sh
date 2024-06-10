#!/bin/bash -l

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --output=./error/simulate_%j.stdout
#SBATCH --error=./error/simulate_%j.stderr
#SBATCH --mail-user=davidhillis@ucsb.edu
#SBATCH --mail-type=ALL
#SBATCH --time=15-00:00:00
#SBATCH --job-name="simulate"
#SBATCH --array=1

# Run in general environment

dir_scripts="/home/davidhillis/anaconda3/envs/general/bin/"

echo "Start 01"
${dir_scripts}Rscript 01_HR_response_to_selection_simulations.R

echo "Start 02"
${dir_scripts}Rscript 02_graph_single_iteration.R

echo "Start 03"
${dir_scripts}Rscript 03_graph_multiple_iterations.R

echo "Start 04"
${dir_scripts}Rscript 04_calculate_power_for_multiple_iterations.R

echo "Start 05"
${dir_scripts}Rscript 05_power_comparisons_between_g22_and_g61.R

echo "Start 06"
${dir_scripts}Rscript 06_power_by_generation.R

echo "Start 07"
${dir_scripts}Rscript 07_power_by_allele_effect.R

echo "Start 08"
${dir_scripts}Rscript 08_selection_differentials.R

echo "Start 09"
${dir_scripts}Rscript 09_heritability.R

