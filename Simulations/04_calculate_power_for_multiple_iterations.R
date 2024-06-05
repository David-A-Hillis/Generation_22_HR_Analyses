# This code was written for: 
# David A Hillis, Liran Yadgary, George M. Weinstock, Fernando Pardo-Manuel
# de Villena, Daniel Pomp, and Theodore Garland Jr. “Large Changes in Detected Selection Signatures After a
# Selection Limit in Mice Bred for Voluntary Wheel-Running Behavior.” PLOS ONE, 2024.

# This calculates the power for the results of 01_HR_response_to_selection_simulations.R and generates a table

data_dir <- "/home/davidhillis/Simulation_trial/Data/"
results_dir <- "/home/davidhillis/Simulation_trial/Results/"
letters <- LETTERS[c(1)] # This defaults to "A" but can include more to separate the results of parallel processing
total_replications <- 10
iter <- c(seq(0,20,5), 22, seq(25,60,5), 61) # These are the generations for which power will be tested
uncon <- 50000
con <- 10000
unconstrained_replications <- 5


gen <- c()
for (i in 1:length(iter)) {gen <- c(gen, paste("gen", iter[i], sep = ""))}
final <- as.data.frame(matrix(NA, ncol = length(iter)+1, nrow = 0))
colnames(final) <- c("constraint", gen)

for (i in letters) {
  for (j in 1:total_replications) {
    setwd(data_dir)
    AF_noconstraint <- NULL
    names <- c()
    
    # Compile data from each of the lines and generations of interest
    for (k in iter) {
      for (p in 1:8) {
        temp <- read.csv(paste("constraint_simulation_line", p, "_2096loci_gen", k, "_run", i, j, ".csv", sep = ""))
        for (m in 1:nrow(temp)) {for (n in 1:ncol(temp)) {
          if (temp[m,n]<0) {temp[m,n] <- 0}
          else {temp[m,n] <- 1}}} 
        temp$af <- rowMeans(temp)
        AF_noconstraint <- cbind(AF_noconstraint, temp[,ncol(temp)]) 
        names <- c(names, paste("L", p, "_g", k, sep = ""))
      }
    }
    
    AF_noconstraint <- as.data.frame(AF_noconstraint)
    colnames(AF_noconstraint) <- names
    pvalues <- as.data.frame(matrix(NA, nrow = nrow(AF_noconstraint), ncol = length(iter)))
    colnames(pvalues) <- gen
    
    # Conduct statistical tests for each locus and generation
    for (k in 1:ncol(pvalues)) {
      for (m in 1:nrow(pvalues)) {
        if (sum(AF_noconstraint[m,c((((k-1)*8)+1):((k*8)-4))]==AF_noconstraint[m,c((((k-1)*8)+5):(k*8))])==4) {pvalues[m,k] <- 0.9999999}
        else if (sum(abs(AF_noconstraint[m,c((((k-1)*8)+1):((k*8)-4))]-AF_noconstraint[m,c((((k-1)*8)+5):(k*8))]))==4 &
                 abs(sum(AF_noconstraint[m,c((((k-1)*8)+1):((k*8)-4))]))==0) {pvalues[m,k] <- 0.00001} 
        else if (sum(abs(AF_noconstraint[m,c((((k-1)*8)+1):((k*8)-4))]-AF_noconstraint[m,c((((k-1)*8)+5):(k*8))]))==4 &
                 abs(sum(AF_noconstraint[m,c((((k-1)*8)+1):((k*8)-4))]))==4) {pvalues[m,k] <- 0.00001} 
        else {test <- t.test(asin(sqrt(AF_noconstraint[m,c((((k-1)*8)+1):((k*8)-4))])), asin(sqrt(AF_noconstraint[m,c((((k-1)*8)+5):(k*8))])))
        pvalues[m,k] <- test$p.value}
      }
    }
    
    if (j<=unconstrained_replications) {filename <- paste("constraint_simulation_2096loci_50AF_unconstrained_run", i, j, sep = "")}
    else if (j>unconstrained_replications) {filename <- paste("constraint_simulation_2096loci_50AF_constraint", con, "_run", i, j, sep = "")}
    
    setwd(results_dir)
    write.csv(AF_noconstraint, paste("simulation_", i, j, "_AF_", filename, ".csv", sep = ""), row.names = FALSE)
    write.csv(pvalues, paste("simulation_", i, j, "_pvalues_", filename, ".csv", sep = ""), row.names = FALSE)
    
    finaltemp <- c()
    for (k in 1:ncol(pvalues)) {finaltemp <- c(finaltemp, sum(pvalues[,k]<=0.05, na.rm = TRUE)/nrow(pvalues))}
    if (j<=unconstrained_replications) {finaltemp <- c("unconstrained", finaltemp)}
    else if (j>unconstrained_replications) {finaltemp <- c("constraint", finaltemp)}
    finaltemp <- as.data.frame(t(as.data.frame(finaltemp)))
    final <- rbind(final, finaltemp)
  }
}

names <- c()
for (i in letters) {for (j in 1:total_replications) {names <- c(names, paste("run", i, j, sep = ""))}}
names <- as.data.frame(names)
final <- cbind(names, final)
colnames(final)[1] <- "Iteration"

write.csv(final, paste("2096_Loci_multi_power", letters[1], ".csv", sep = ""), row.names = FALSE)
