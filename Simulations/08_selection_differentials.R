# This code was written for: 
# David A Hillis, Liran Yadgary, George M. Weinstock, Fernando Pardo-Manuel
# de Villena, Daniel Pomp, and Theodore Garland Jr. “Large Changes in Detected Selection Signatures After a
# Selection Limit in Mice Bred for Voluntary Wheel-Running Behavior.” PLOS ONE, 2024.

# This calculates selection differentials for the results of 01_HR_response_to_selection_simulations.R 

data_dir <- "/home/davidhillis/Simulation_trial/Data/"
results_dir <- "/home/davidhillis/Simulation_trial/Results/"

indLines <- 1 # set to 0 for means, set to 1 for individual lines
letters <- LETTERS[c(1)]
unconstrained_replications <- 5
total_replications <- 10

final <- as.data.frame(matrix(NA, nrow = 61, ncol = 0))

# HR selection differential unconstrained
setwd(data_dir)
for (m in letters) {
  for (n in 1:unconstrained_replications) {
    breed <- read.csv(paste("constraint_simulation_2096loci_50AF_v77_constraint50000_run", m, n, ".csv", sep = ""), 
                      stringsAsFactors = FALSE)
    breed <- cbind(breed[,1], as.data.frame(rep(c(0:61), 4)), breed[,c(2:ncol(breed))])
    data <- read.csv(paste("wholeRunning_HR_2096loci_50AF_v77_constraint50000_run", m, n, ".csv", sep = ""), 
                     stringsAsFactors = FALSE)
    data <- cbind(data[,1], as.data.frame(rep(c(1:61), 4)), data[,c(2:ncol(data))])
    
    final_SelDif <- as.data.frame(matrix(NA, nrow = 61, ncol = 80))
    
    for (i in 5:8) {
      data_temp <- data[which(data[,1]==i),]
      data_temp <- data_temp[,-c(1,2)]
      for (j in 1:61) {
        SD <- sd(c(as.numeric(data[which(data[,2]==j & data[,1]==i),-c(1,2)])))
        breed_temp <- breed[which(breed[,2]==j & breed[,1]==i),-c(1,2)]
        for (k in 1:20) {
          group <- as.numeric(data_temp[j, ((((k-1)*5)+1):(((k-1)*5)+5))])
          final_SelDif[j,(((i-5)*20)+k)] <-  (breed_temp[k]-mean(group))/sd(group) # SD
        }
      }
    }
    if (indLines==0) {final <- cbind(final, rowMeans(final_SelDif))}
    else if (indLines==1) {for (i in 1:4) {final <- cbind(final, rowMeans(final_SelDif[,c((((i-1)*20)+1):(((i-1)*20)+20))]))}}
  }
}

setwd(results_dir)
if (indLines==0) {jpeg(filename = paste("Selection_differential_unconstrained_allsims_SDgroup_v100.jpeg", sep = ""), width = 4600, height = 2300, res = 500)}
if (indLines==1) {jpeg(filename = paste("Selection_differential_unconstrained_allsims_SDgroup_v100_indlines.jpeg", sep = ""), width = 4600, height = 2300, res = 500)}
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1:61), final[,1], type = "l", lwd = 2, axes = FALSE, ylim = c(0,1.4), cex.main = 1.2, cex.lab = 1,
     main = "", xlab = "Generation", ylab = "Selection Differential")
for (i in 2:ncol(final)) {lines(c(1:61), final[,i], lwd = 2)}
lines(c(1:61), rowMeans(final, na.rm = TRUE), lwd = 2, col = "red")
axis(side = 1, seq(1, 61, 2), cex.axis=1)
axis(side = 2, seq(0,1.4, 0.1), cex.axis=1)
dev.off()


# HR selection differential constrained
# With default settings, the graph contains dips for seasonal change
letters <- LETTERS[c(1)]
setwd(data_dir)
final <- as.data.frame(matrix(NA, nrow = 61, ncol = 0))

for (m in letters) {
  for (n in c((unconstrained_replications+1):total_replications)) {
    breed <- read.csv(paste("constraint_simulation_2096loci_50AF_v77_constraint10000_run", m, n, ".csv", sep = ""), 
                      stringsAsFactors = FALSE)
    breed <- cbind(breed[,1], as.data.frame(rep(c(0:61), 4)), breed[,c(2:ncol(breed))])
    data <- read.csv(paste("wholeRunning_HR_2096loci_50AF_v77_constraint10000_run", m, n, ".csv", sep = ""), 
                     stringsAsFactors = FALSE)
    data <- cbind(data[,1], as.data.frame(rep(c(1:61), 4)), data[,c(2:ncol(data))])
    
    final_SelDif <- as.data.frame(matrix(NA, nrow = 61, ncol = 80))
    for (i in 1:nrow(breed)) {for (k in 1:ncol(breed)) {if (breed[i,k]>=16000) {breed[i,k] <- 16000}}}
    for (i in 1:nrow(data)) {for (k in 1:ncol(data)) {if (data[i,k]>=16000) {data[i,k] <- 16000}}}
    
    for (i in 5:8) {
      data_temp <- data[which(data[,1]==i),]
      data_temp <- data_temp[,-c(1,2)]
      for (j in 1:61) {
        #SD <- sd(c(as.numeric(data[which(data[,2]==j & data[,1]==i),-c(1,2)])))
        breed_temp <- breed[which(breed[,2]==j & breed[,1]==i),-c(1,2)]
        for (k in 1:20) {
          group <- as.numeric(data_temp[j, ((((k-1)*5)+1):(((k-1)*5)+5))])
          if (sd(group)==0) {final_SelDif[j,(((i-5)*20)+k)] <- 0}
          else {final_SelDif[j,(((i-5)*20)+k)] <-  (breed_temp[k]-mean(group))/sd(group)} # SD
        }
      }
    }
    if (indLines==0) {final <- cbind(final, rowMeans(final_SelDif))}
    else if (indLines==1) {for (i in 1:4) {final <- cbind(final, rowMeans(final_SelDif[,c((((i-1)*20)+1):(((i-1)*20)+20))]))}}
  }
}

setwd(results_dir)
if (indLines==0) {jpeg(filename = paste("Selection_differential_constrained16000_allsims_SDgroup_set0_v100.jpeg", sep = ""), width = 4600, height = 2300, res = 500)}
if (indLines==1) {jpeg(filename = paste("Selection_differential_constrained16000_allsims_SDgroup_set0_v100_indlines.jpeg", sep = ""), width = 4600, height = 2300, res = 500)}
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1:61), final[,1], type = "l", lwd = 2, axes = FALSE, ylim = c(0,1.4), cex.main = 1.2, cex.lab = 1,
     main = "", xlab = "Generation", ylab = "Selection Differential")
for (i in 2:ncol(final)) {lines(c(1:61), final[,i], lwd = 2)}
lines(c(1:61), rowMeans(final, na.rm = TRUE), lwd = 2, col = "red")
axis(side = 1, seq(1, 61, 2), cex.axis=1)
axis(side = 2, seq(0,1.4, 0.1), cex.axis=1)
dev.off()

# HR selection differential constrained
# This does not include actual running levels but instead the amount each mouse were to run 
# if the constraint were suddenly removed
letters <- LETTERS[c(1)]
setwd(data_dir)
final <- as.data.frame(matrix(NA, nrow = 61, ncol = 0))

for (m in letters) {
  for (n in c((unconstrained_replications+1):total_replications)) {
    breed <- read.csv(paste("constraint_simulation_2096loci_50AF_v77_constraint10000_run", m, n, ".csv", sep = ""), 
                      stringsAsFactors = FALSE)
    breed <- cbind(breed[,1], as.data.frame(rep(c(0:61), 4)), breed[,c(2:ncol(breed))])
    data <- read.csv(paste("wholeRunning_HR_2096loci_50AF_v77_constraint10000_run", m, n, ".csv", sep = ""), 
                     stringsAsFactors = FALSE)
    data <- cbind(data[,1], as.data.frame(rep(c(1:61), 4)), data[,c(2:ncol(data))])
    
    final_SelDif <- as.data.frame(matrix(NA, nrow = 61, ncol = 80))
    
    for (i in 5:8) {
      data_temp <- data[which(data[,1]==i),]
      data_temp <- data_temp[,-c(1,2)]
      for (j in 1:61) {
        # SD <- sd(c(as.numeric(data[which(data[,2]==j & data[,1]==i),-c(1,2)])))
        breed_temp <- breed[which(breed[,2]==j & breed[,1]==i),-c(1,2)]
        for (k in 1:20) {
          group <- as.numeric(data_temp[j, ((((k-1)*5)+1):(((k-1)*5)+5))])
          if (sd(group)==0) {final_SelDif[j,(((i-5)*20)+k)] <- 0}
          else {final_SelDif[j,(((i-5)*20)+k)] <-  (breed_temp[k]-mean(group))/sd(group)} # SD
        }
      }
    }
    if (indLines==0) {final <- cbind(final, rowMeans(final_SelDif))}
    else if (indLines==1) {for (i in 1:4) {final <- cbind(final, rowMeans(final_SelDif[,c((((i-1)*20)+1):(((i-1)*20)+20))]))}}
  }
}

setwd(results_dir)
if (indLines==0) {jpeg(filename = paste("Selection_differential_constrained10000_allsims_notreal_SDgroup_set0_v100.jpeg", sep = ""), width = 4600, height = 2300, res = 500)}
if (indLines==1) {jpeg(filename = paste("Selection_differential_constrained10000_allsims_notreal_SDgroup_set0_v100_indlines.jpeg", sep = ""), width = 4600, height = 2300, res = 500)}
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1:61), final[,1], type = "l", lwd = 2, axes = FALSE, ylim = c(0,1.4), cex.main = 1.2, cex.lab = 1,
     main = "", xlab = "Generation", ylab = "Selection Differential")
for (i in 2:ncol(final)) {lines(c(1:61), final[,i], lwd = 2)}
lines(c(1:61), rowMeans(final, na.rm = TRUE), lwd = 2, col = "red")
axis(side = 1, seq(1, 61, 2), cex.axis=1)
axis(side = 2, seq(0,1.4, 0.1), cex.axis=1)
dev.off()


# Control selection differential 
indLines <- 0 # set to 0 for means, set to 1 for individual lines
letters <- LETTERS[c(1)] 
setwd(data_dir)
final <- as.data.frame(matrix(NA, nrow = 61, ncol = 0))

for (m in letters) {
  for (n in c(1:total_replications)) {
    if (n<=unconstrained_replications) {con <- 50000}
    else {con <- 10000}
    breed <- read.csv(paste("constraint_simulation_2096loci_50AF_v77_constraint", con, "_run", m, n, ".csv", sep = ""), 
                      stringsAsFactors = FALSE)
    breed <- cbind(breed[,1], as.data.frame(rep(c(0:61), 4)), breed[,c(2:ncol(breed))])
    data <- read.csv(paste("wholeRunning_C_2096loci_50AF_v77_constraint", con, "_run", m, n, ".csv", sep = ""), 
                     stringsAsFactors = FALSE)
    data <- cbind(data[,1], as.data.frame(rep(c(1:61), 4)), data[,c(2:ncol(data))])
    
    final_SelDif <- as.data.frame(matrix(NA, nrow = 61, ncol = 80))
    
    for (i in 1:4) {
      data_temp <- data[which(data[,1]==i),]
      data_temp <- data_temp[,-c(1,2)]
      for (j in 1:61) {
        # SD <- sd(c(as.numeric(data[which(data[,2]==j & data[,1]==i),-c(1,2)])))
        breed_temp <- breed[which(breed[,2]==j & breed[,1]==i),-c(1,2)]
        for (k in 1:20) {
          group <- as.numeric(data_temp[j, ((((k-1)*2)+1):(((k-1)*2)+2))])
          if (sd(group)==0) {final_SelDif[j,(((i-1)*20)+k)] <- 0}
          else {final_SelDif[j,(((i-1)*20)+k)] <-  (breed_temp[k]-mean(group))/sd(group)} # SD
        }
      }
    }
    if (indLines==0) {final <- cbind(final, rowMeans(final_SelDif))}
    else if (indLines==1) {for (i in 1:4) {final <- cbind(final, rowMeans(final_SelDif[,c((((i-1)*20)+1):(((i-1)*20)+20))]))}}
  }
}

setwd(results_dir)
if (indLines==0) {jpeg(filename = paste("Selection_differential_control_v100_half.jpeg", sep = ""), width = 4600, height = 2300, res = 500)}
if (indLines==1) {jpeg(filename = paste("Selection_differential_control_v100_indlines.jpeg", sep = ""), width = 4600, height = 2300, res = 500)}
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1:61), final[,1], type = "l", lwd = 2, axes = FALSE, ylim = c(-0.7,0.7), cex.main = 1.2, cex.lab = 1,
     main = "", xlab = "Generation", ylab = "Selection Differential")
for (i in 2:ncol(final)) {lines(c(1:61), final[,i], lwd = 2)}
lines(c(1:61), rowMeans(final, na.rm = TRUE), lwd = 2, col = "red")
axis(side = 1, seq(1, 61, 2), cex.axis=1)
axis(side = 2, seq(-0.7,0.7, 0.1), labels = seq(-0.7,0.7, 0.1), cex.axis=1)
dev.off()



