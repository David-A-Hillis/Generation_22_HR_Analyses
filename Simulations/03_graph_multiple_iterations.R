# This code was written for: 
# David A Hillis, Liran Yadgary, George M. Weinstock, Fernando Pardo-Manuel
# de Villena, Daniel Pomp, and Theodore Garland Jr. “Large Changes in Detected Selection Signatures After a
# Selection Limit in Mice Bred for Voluntary Wheel-Running Behavior.” PLOS ONE, 2024.

# This graphs the multiple iterations of results of the simulations produced by 01_HR_response_to_selection_simulations.R


data_dir <- "/home/davidhillis/Simulation_trial/Data/"
results_dir <- "/home/davidhillis/Simulation_trial/Results/"

total_replications <- 10 # Number of iterations for each letter
gen <- 61 # Number of generations
uncon <- 50000 # Unconstrained limit
con <- 10000 # Constrained limit
unconstrained_replications <- 5 # How many of the total replications were unconstrained

for (j in 1:total_replications) {
  setwd(data_dir)
  if (j<=unconstrained_replications) {constraint <- uncon} # 
  else if (j>unconstrained_replications) {constraint <- con}
  runningC <- read.csv(paste("wholeRunning_C_2096loci_50AF_v77_constraint", constraint,"_runA", j, ".csv", sep = ""), stringsAsFactors = FALSE)
  runC <- runningC[,c(1, seq(2,38,4), seq(3,39,4))] # one sex
  runningHR <- read.csv(paste("wholeRunning_HR_2096loci_50AF_v77_constraint", constraint,"_runA", j, ".csv", sep = "")
                        , stringsAsFactors = FALSE)
  runHR <- runningHR[,c(1, seq(2,92,10), seq(3,93,10), seq(4,94,10), seq(5,95,10), seq(6,96,10))] # one sex
  
  runC$mean <- rowMeans(runC[,c(2:ncol(runC))])
  runHR$mean <- rowMeans(runHR[,c(2:ncol(runHR))])
  running <- rbind(runC[,c(1,ncol(runC))], runHR[,c(1,ncol(runHR))])
  
  setwd(results_dir)
  if (j<=unconstrained_replications) {jpeg(filename = paste("unconstrainedA", j, ".jpg", sep = ""), width = 4600, height = 2300, res = 500)}
  else if (j>unconstrained_replications) {jpeg(filename = paste("constrainedA", j, "_real.jpg", sep = ""), width = 4600, height = 2300, res = 500)}
  par(mar = c(5, 6, 4, 2) + 0.1) #par(mar = c(10, 10, 4, 2) + 0.1, mgp = c(6,1,0))
  plot(seq(1, gen, 1), running[which(running[,1]==1),2], type = "l", col = "dodgerblue1", lwd = 2, ylim = c(0,30000), 
       xlab = "Generation", ylab = "Revolutions per Day")
  for (k in 1:15) {abline(h=k*2000, col = "gray73", lwd = 2)}
  lines(seq(1, gen, 1), running[which(running[,1]==2),2], type = "l", col = "dodgerblue2", lwd = 2)
  lines(seq(1, gen, 1), running[which(running[,1]==3),2], type = "l", col = "dodgerblue3", lwd = 2)
  lines(seq(1, gen, 1), running[which(running[,1]==4),2], type = "l", col = "dodgerblue4", lwd = 2)
  lines(seq(1, gen, 1), running[which(running[,1]==5),2], type = "l", col = "firebrick1", lwd = 2)
  lines(seq(1, gen, 1), running[which(running[,1]==6),2], type = "l", col = "firebrick2", lwd = 2)
  lines(seq(1, gen, 1), running[which(running[,1]==7),2], type = "l", col = "firebrick3", lwd = 2)
  lines(seq(1, gen, 1), running[which(running[,1]==8),2], type = "l", col = "firebrick4", lwd = 2)
  dev.off()
  
}
