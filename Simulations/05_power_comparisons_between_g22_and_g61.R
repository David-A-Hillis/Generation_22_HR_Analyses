# This code was written for: 
# David A Hillis, Liran Yadgary, George M. Weinstock, Fernando Pardo-Manuel
# de Villena, Daniel Pomp, and Theodore Garland Jr. “Large Changes in Detected Selection Signatures After a
# Selection Limit in Mice Bred for Voluntary Wheel-Running Behavior.” PLOS ONE, 2024.

# This compares the power of gen22 and gen61 for the results of 01_HR_response_to_selection_simulations.R 
# Note these results are not written to the Results folder, but are printed within R 

results_dir <- "/home/davidhillis/Simulation_trial/Results/"
letters <- LETTERS[c(1)] # This defaults to "A" but can include more to separate the results of parallel processing
total_replications <- 10
unconstrained_replications <- 5

final_con <- c()
final_uncon <- c()

setwd(results_dir)

# Combine data and perform arcsine transformation to prep for statistical tests
for (i in letters) {
  for (j in 1:total_replications) {
    if (j<=unconstrained_replications) {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_unconstrained_run", i, j, ".csv", sep = ""))
      final_uncon <- c(final_uncon, cor(asin(sqrt(temp[,6])), asin(sqrt(temp[,15]))))}
    else if (j>unconstrained_replications & i!="G" & i!="P") {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_constraint10000_run", i, j, ".csv", sep = ""))
      final_con <- c(final_con, cor(asin(sqrt(temp[,6])), asin(sqrt(temp[,15]))))}
    else if (i=="G" | i=="P") {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_constraint10000_run", i, j, ".csv", sep = ""))
      final_con <- c(final_con, cor(asin(sqrt(temp[,6])), asin(sqrt(temp[,15]))))}
  }
}

# Comparing the correlation between generations 22 and 61 for all 2096 loci
t.test(final_con, final_uncon)

final_con <- c()
final_uncon <- c()

for (i in letters) {
  for (j in 1:total_replications) {
    if (j<=unconstrained_replications) {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_unconstrained_run", i, j, ".csv", sep = ""))
      sig <- temp[,6]<=0.05
      sig2 <- temp[,15]<=0.05
      final_uncon <- c(final_uncon, sum(sig==1 & sig==sig2))}
    else if (j>unconstrained_replications & i!="G" & i!="P") {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_constraint10000_run", i, j, ".csv", sep = ""))
      sig <- temp[,6]<=0.05
      sig2 <- temp[,15]<=0.05
      final_con <- c(final_con, sum(sig==1 & sig==sig2))}
    else if (i=="G" | i=="P") {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_constraint10000_run", i, j, ".csv", sep = ""))
      sig <- temp[,6]<=0.05
      sig2 <- temp[,15]<=0.05
      final_con <- c(final_con, sum(sig==1 & sig==sig2))}
  }
}

# Comparing the number of agreed positives between generations 20 and 60
t.test(final_con, final_uncon)

final_con <- c()
final_uncon <- c()

for (i in letters) {
  for (j in 1:total_replications) {
    if (j<=unconstrained_replications) {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_unconstrained_run", i, j, ".csv", sep = ""))
      sig <- temp[,6]<=0.05
      sig2 <- temp[,15]<=0.05
      final_uncon <- c(final_uncon, (sum(sig==1 & sig==sig2)/sum(sig==1)))}
    else if (j>unconstrained_replications & i!="G" & i!="P") {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_constraint10000_run", i, j, ".csv", sep = ""))
      sig <- temp[,6]<=0.05
      sig2 <- temp[,15]<=0.05
      final_con <- c(final_con, (sum(sig==1 & sig==sig2)/sum(sig==1)))}
    else if (i=="G" | i=="P") {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_constraint10000_run", i, j, ".csv", sep = ""))
      sig <- temp[,6]<=0.05
      sig2 <- temp[,15]<=0.05
      final_con <- c(final_con, (sum(sig==1 & sig==sig2)/sum(sig==1)))}
  }
}

# Comparing the proportion of agreed positives between generations 20 and 60
t.test(final_con, final_uncon)

