# This code was written for: 
# David A Hillis, Liran Yadgary, George M. Weinstock, Fernando Pardo-Manuel
# de Villena, Daniel Pomp, and Theodore Garland Jr. “Large Changes in Detected Selection Signatures After a
# Selection Limit in Mice Bred for Voluntary Wheel-Running Behavior.” PLOS ONE, 2024.

# This calculates power by allele effect size for the results of 01_HR_response_to_selection_simulations.R 

results_dir <- "/home/davidhillis/Simulation_trial/Results/"
letters <- LETTERS[c(1)]
unconstrained_replications <- 5
total_replications <- 10

# Allelic weights
gweight <- c(rep(0.4, 720), rep(0.8, 480), rep(1.6, 312), rep(3.2, 216), rep(6.4, 144), rep(12.8, 96), rep(25.6, 60), rep(51.2, 36), rep(102.4, 24), rep(204.8, 8)) 

# Generate output dataframes
finalcon22 <- as.data.frame(matrix(NA, ncol = 0, nrow = 10))
finaluncon22 <- as.data.frame(matrix(NA, ncol = 0, nrow = 10))
finalcon61 <- as.data.frame(matrix(NA, ncol = 0, nrow = 10))
finaluncon61 <- as.data.frame(matrix(NA, ncol = 0, nrow = 10))

# Aggregate data
setwd(results_dir)
for (i in letters) {
  for (j in 1:total_replications) {
    if (j<=unconstrained_replications) {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_unconstrained_run", i, j, ".csv", sep = ""))
      temp <- cbind(as.data.frame(gweight), temp)
      sig <- c()
      for (k in unique(gweight)) {sig <- c(sig, nrow(temp[which(temp[,7]<=0.05 & temp[,1]==k),]))}
      finaluncon22 <- cbind(finaluncon22, as.data.frame(sig))
      
      sig <- c()
      for (k in unique(gweight)) {sig <- c(sig, nrow(temp[which(temp[,16]<=0.05 & temp[,1]==k),]))}
      finaluncon61 <- cbind(finaluncon61, as.data.frame(sig))}
    
    else if (j>unconstrained_replications) {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_constraint10000_run", i, j, ".csv", sep = ""))
      temp <- cbind(as.data.frame(gweight), temp)
      sig <- c()
      for (k in unique(gweight)) {sig <- c(sig, nrow(temp[which(temp[,7]<=0.05 & temp[,1]==k),]))}
      finalcon22 <- cbind(finalcon22, as.data.frame(sig))
      
      sig <- c()
      for (k in unique(gweight)) {sig <- c(sig, nrow(temp[which(temp[,16]<=0.05 & temp[,1]==k),]))}
      finalcon61 <- cbind(finalcon61, as.data.frame(sig))}
  }
}

setwd(results_dir)
final <- as.data.frame(matrix(NA, nrow = 10, ncol = 13))
colnames(final) <- c("Affect_size", "unconstrained_g22_mean", "unconstrained_g22_power", "constrained_g22_mean", "constrained_g22_power", 
                     "unconstrained_g61_mean", "unconstrained_g61_power", "constrained_g61_mean", "constrained_g61_power", 
                     "pvalue_g22", "pvalue_g61", "pvalue_unconstrained", "pvalue_constrained")
                     
# Conduct statistical tests for differences between generations 22 and 61
n <- 1
for (i in unique(gweight)) {
  final[n,1] <- i
  final[n,2] <- mean(as.numeric(finaluncon22[n,]))
  final[n,3] <- mean(as.numeric(finaluncon22[n,]))/nrow(temp[which(temp[,1]==i),])
  final[n,4] <- mean(as.numeric(finalcon22[n,]))
  final[n,5] <- mean(as.numeric(finalcon22[n,]))/nrow(temp[which(temp[,1]==i),])
  final[n,6] <- mean(as.numeric(finaluncon61[n,]))
  final[n,7] <- mean(as.numeric(finaluncon61[n,]))/nrow(temp[which(temp[,1]==i),])
  final[n,8] <- mean(as.numeric(finalcon61[n,]))
  final[n,9] <- mean(as.numeric(finalcon61[n,]))/nrow(temp[which(temp[,1]==i),])
  test <- t.test(as.numeric(finaluncon22[n,]), as.numeric(finalcon22[n,]))
  final[n,10] <- test$p.value
  test <- t.test(as.numeric(finaluncon61[n,]), as.numeric(finalcon61[n,]))
  final[n,11] <- test$p.value
  test <- t.test(as.numeric(finaluncon22[n,]), as.numeric(finaluncon61[n,]))
  final[n,12] <- test$p.value
  test <- t.test(as.numeric(finalcon22[n,]), as.numeric(finalcon61[n,]))
  final[n,13] <- test$p.value
  n <- n+1
}

write.csv(final, "Power_comparisons_by_affect_sizes.csv", row.names = FALSE)

# For generation 0

finalcon0 <- as.data.frame(matrix(NA, ncol = 0, nrow = 10))
finaluncon0 <- as.data.frame(matrix(NA, ncol = 0, nrow = 10))

# The code below is a repitition of what is above, except for generation 0, for comparison  
setwd(results_dir)
for (i in letters) {
  for (j in 1:total_replications) {
    if (j<=unconstrained_replications) {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_unconstrained_run", i, j, ".csv", sep = ""))
      temp <- cbind(as.data.frame(gweight), temp)
      sig <- c()
      for (k in unique(gweight)) {sig <- c(sig, nrow(temp[which(temp[,2]<=0.05 & temp[,1]==k),]))}
      finaluncon0 <- cbind(finaluncon0, as.data.frame(sig))}
    
    else if (j>unconstrained_replications) {
      temp <- read.csv(paste("simulation_", i, j, "_pvalues_constraint_simulation_2096loci_50AF_constraint10000_run", i, j, ".csv", sep = ""))
      temp <- cbind(as.data.frame(gweight), temp)
      sig <- c()
      for (k in unique(gweight)) {sig <- c(sig, nrow(temp[which(temp[,2]<=0.05 & temp[,1]==k),]))}
      finalcon0 <- cbind(finalcon0, as.data.frame(sig))}
  }
}

final <- as.data.frame(matrix(NA, nrow = 10, ncol = 6))
colnames(final) <- c("Affect_size", "unconstrained_g0_mean", "unconstrained_g0_power", "constrained_g0_mean", "constrained_g0_power", "pvalue_g0")
n <- 1
for (i in unique(gweight)) {
  final[n,1] <- i
  final[n,2] <- mean(as.numeric(finaluncon0[n,]))
  final[n,3] <- mean(as.numeric(finaluncon0[n,]))/nrow(temp[which(temp[,1]==i),])
  final[n,4] <- mean(as.numeric(finalcon0[n,]))
  final[n,5] <- mean(as.numeric(finalcon0[n,]))/nrow(temp[which(temp[,1]==i),])
  test <- t.test(as.numeric(finaluncon0[n,]), as.numeric(finalcon0[n,]))
  final[n,6] <- test$p.value
  n <- n+1
}

setwd(results_dir)
write.csv(final, "Power_comparisons_by_affect_sizes_gen0.csv", row.names = FALSE)

