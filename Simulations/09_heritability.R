# This code was written for: 
# David A Hillis, Liran Yadgary, George M. Weinstock, Fernando Pardo-Manuel
# de Villena, Daniel Pomp, and Theodore Garland Jr. “Large Changes in Detected Selection Signatures After a
# Selection Limit in Mice Bred for Voluntary Wheel-Running Behavior.” PLOS ONE, 2024.

# This calculates heritability for the results of 01_HR_response_to_selection_simulations.R 

data_dir <- "/home/davidhillis/Simulation_trial/Data/"
results_dir <- "/home/davidhillis/Simulation_trial/Results/"

gen <- 61 # Number of generations in the simulation
letters <- LETTERS[c(1)] 
unconstrained_replications <- 5
total_replications <- 10
clines <- 4 # Number of control lines
cmice <- 2 # Number of pups per family in control lines
hrlines <- 4 # Number of HR lines
hrmice <- 5 # Number of pups per family in HR lines

# Create output dataframes
final_hetC <- as.data.frame(matrix(NA, ncol = 0, nrow = gen))
final_hetHRcon <- as.data.frame(matrix(NA, ncol = 0, nrow = gen))
final_hetHRuncon <- as.data.frame(matrix(NA, ncol = 0, nrow = gen))


for (m in letters) {
  for (n in 1:total_replications) {
    setwd(data_dir)
    if (n>unconstrained_replications) {
      running <- read.csv(paste("constraint_simulation_2096loci_50AF_v77_constraint10000_run", m, n, ".csv", sep = ""), stringsAsFactors = FALSE)
      for (i in 2:21) {for (j in 248:496) {if (running[j,i]>16000) {running[j,i]<-10000}}}}
    else if (n<=unconstrained_replications) {
      running <- read.csv(paste("constraint_simulation_2096loci_50AF_v77_constraint50000_run", m, n, ".csv", sep = ""), stringsAsFactors = FALSE)}
    
    final_het <- as.data.frame(matrix(NA, ncol = (clines+hrlines), nrow = 0))
    
    cnames <- c()
    for (i in 1:clines) {cnames <- c(cnames, paste("C", i))}
    for (i in (clines+1):(clines+hrlines)) {cnames <- c(cnames, paste("HR", i))}
    colnames(final_het) <- cnames
    
    # control heritability
    if (n>unconstrained_replications) {hr_run <- read.csv(paste("wholeRunning_C_2096loci_50AF_v77_constraint10000_run", m, n, ".csv", sep = ""), stringsAsFactors = FALSE)}
    else if (n<=unconstrained_replications) {hr_run <- read.csv(paste("wholeRunning_C_2096loci_50AF_v77_constraint50000_run", m, n, ".csv", sep = ""), stringsAsFactors = FALSE)}
    hr_run <- cbind(hr_run[,1], as.data.frame(rep(c(1:gen), clines)), hr_run[,c(2:ncol(hr_run))])
    for (i in 1:gen) {
      for (j in 1:clines) {
        pairs <- read.csv(paste("breedOrder_line", j, "_2096loci_gen", i, "_run", m, n, ".csv", sep = ""), stringsAsFactors = FALSE)
        data_temp <- as.data.frame(matrix(NA, ncol = 2, nrow = 10))
        colnames(data_temp) <- c("midparent", "midoffspring")
        running_temp <- running[which(running[,1]==j),]
        running_temp <- running_temp[(i:(i+1)),-1]
        hr_temp <- hr_run[which(hr_run[,1]==j & hr_run[,2]==i),-c(1,2)]
        
        for (k in 1:10) {
          data_temp[k,1] <- mean(c(running_temp[1,(((k-1)*2)+1)], running_temp[1,(((pairs[k,2]-1)*2)+2)]))
          data_temp[k,2] <- mean(as.numeric(hr_temp[1,c((((k-1)*(cmice*2))+1):(((k-1)*(cmice*2))+(cmice*2)))]))
        }
        lm <- lm(data_temp[,2] ~ data_temp[,1])
        final_het[i,j] <- lm$coefficients[2]
      }
    }
    
    # HR heritability
    if (n>unconstrained_replications) {
      hr_run <- read.csv(paste("wholeRunning_HR_2096loci_50AF_v77_constraint10000_run", m, n, ".csv", sep = ""), stringsAsFactors = FALSE)
      for (i in 2:101) {for (j in 1:244) {if (hr_run[j,i]>16000) {hr_run[j,i]<-16000}}}}
    else if (n<=unconstrained_replications) {hr_run <- read.csv(paste("wholeRunning_HR_2096loci_50AF_v77_constraint50000_run", m, n, ".csv", sep = ""), stringsAsFactors = FALSE)}
    
    hr_run <- cbind(hr_run[,1], as.data.frame(rep(c(1:gen), hrlines)), hr_run[,c(2:ncol(hr_run))])
    
    for (i in 1:gen) {
      for (j in (clines+1):(clines+hrlines)) {
        pairs <- read.csv(paste("breedOrder_line", j, "_2096loci_gen", i, "_run", m, n, ".csv", sep = ""), stringsAsFactors = FALSE)
        data_temp <- as.data.frame(matrix(NA, ncol = 2, nrow = 10))
        colnames(data_temp) <- c("midparent", "midoffspring")
        running_temp <- running[which(running[,1]==j),]
        running_temp <- running_temp[(i:(i+1)),-1]
        hr_temp <- hr_run[which(hr_run[,1]==j & hr_run[,2]==i),-c(1,2)]
        
        for (k in 1:total_replications) {
          data_temp[k,1] <- mean(c(running_temp[1,(((k-1)*2)+1)], running_temp[1,(((pairs[k,2]-1)*2)+2)]), na.rm = TRUE)
          data_temp[k,2] <- mean(as.numeric(hr_temp[1,c((((k-1)*(hrmice*2))+1):(((k-1)*(hrmice*2))+(hrmice*2)))]), na.rm = TRUE)
        }
        lm <- lm(data_temp[,2] ~ data_temp[,1])
        final_het[i,j] <- lm$coefficients[2]
      }
    }
    final_het$mean_C <- rowMeans(final_het[,c(1:clines)], na.rm = TRUE)
    final_het$mean_HR <- rowMeans(final_het[,c((clines+1):(clines+hrlines))], na.rm = TRUE)
    
    setwd(results_dir)
    write.csv(final_het, paste("heritability_all_gens_sim", m, n, ".csv", sep=""), row.names = FALSE)
    
    final_hetC <- cbind(final_hetC, final_het[,(ncol(final_het)-1)])
    if (n>unconstrained_replications) {final_hetHRcon <- cbind(final_hetHRcon, final_het[,(ncol(final_het))])}
    if (n<=unconstrained_replications) {final_hetHRuncon <- cbind(final_hetHRuncon, final_het[,(ncol(final_het))])}
  }
}

#hetmin <- min(final_hetC, final_hetHRcon, final_hetHRuncon)
#hetmax <- max(final_hetC, final_hetHRcon, final_hetHRuncon)

setwd(results_dir)
jpeg(filename = paste("Heritability_constrained10000_allsims.jpeg", sep = ""), width = 4600, height = 2300, res = 500)
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1:gen), final_hetHRcon[,1], type = "l", lwd = 2, axes = FALSE, ylim = c(-0.1,1.1), cex.main = 1.2, cex.lab = 1,
     main = "", xlab = "Generation", ylab = "Heritability")
for (i in 2:ncol(final_hetHRcon)) {lines(c(1:61), final_hetHRcon[,i], lwd = 2)}
lines(c(1:gen), rowMeans(final_hetHRcon, na.rm = TRUE), lwd = 2, col = "red")
axis(side = 1, seq(1, 61, 2), cex.axis=1)
axis(side = 2, seq(-0.1,1.1, 0.1), cex.axis=1)
dev.off()

setwd(results_dir)
jpeg(filename = paste("Heritability_unconstrained50000_allsims.jpeg", sep = ""), width = 4600, height = 2300, res = 500)
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1:gen), final_hetHRuncon[,1], type = "l", lwd = 2, axes = FALSE, ylim = c(-0.1,1.1), cex.main = 1.2, cex.lab = 1,
     main = "", xlab = "Generation", ylab = "Heritability")
for (i in 2:ncol(final_hetHRuncon)) {lines(c(1:61), final_hetHRuncon[,i], lwd = 2)}
lines(c(1:gen), rowMeans(final_hetHRuncon, na.rm = TRUE), lwd = 2, col = "red")
axis(side = 1, seq(1, 61, 2), cex.axis=1)
axis(side = 2, seq(-0.1,1.1, 0.1), cex.axis=1)
dev.off()

# Only the first half of C results are used 
# This is to keep the sample size matched with previous 2 graphs
setwd(results_dir)
jpeg(filename = paste("Heritability_control_allsims.jpeg", sep = ""), width = 4600, height = 2300, res = 500)
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1:gen), final_hetC[,1], type = "l", lwd = 2, axes = FALSE, ylim = c(-0.1,1.1), cex.main = 1.2, cex.lab = 1,
     main = "", xlab = "Generation", ylab = "Heritability")
for (i in 2:ceiling(ncol(final_hetC)/2)) {lines(c(1:61), final_hetC[,i], lwd = 2)}
lines(c(1:gen), rowMeans(final_hetC[,c(1:ncol(final_hetHRuncon))], na.rm = TRUE), lwd = 2, col = "red")
axis(side = 1, seq(1, 61, 2), cex.axis=1)
axis(side = 2, seq(-0.1,1.1, 0.1), cex.axis=1)
dev.off()



