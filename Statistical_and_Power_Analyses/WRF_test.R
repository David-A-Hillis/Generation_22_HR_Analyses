# This code was written for: 
# David A Hillis, Liran Yadgary, George M. Weinstock, Fernando Pardo-Manuel
# de Villena, Daniel Pomp, and Theodore Garland Jr. “Large Changes in Detected Selection Signatures After a
# Selection Limit in Mice Bred for Voluntary Wheel-Running Behavior.” PLOS ONE, 2024.

# This code is meant to conduct genomic differentiation analyses between 4 High-Runner and 4 Control lines
# The statistical method used is the windowed, regularized F-test described in the publication.
# This is conducted on both generation 22 and generation 61 allele frequencies.

generation <- 61 # make equal to 22 or 61 depending on generation to be analyzed
dir_in <- "/home/davidhillis/Simulation_trial/WRF/"
dir_out <- "/home/davidhillis/Simulation_trial/WRF/"

setwd(dir_in)
data_raw <- read.csv(paste("gen", generation, "_all_AF.csv", sep = ""), stringsAsFactors = FALSE)

chr <- unique(data_raw[,1]) # will probably just eliminate everything that isn't 1-19
data_raw <- data_raw[order(data_raw[,1], data_raw[,2]),]

temp <- data_raw[,c(1:2)]
for (i in 1:8) {temp[,i+2] <- asin(sqrt(data_raw[,i+2]))} # applies arcsine transform for allelic data
colnames(temp) <- colnames(data_raw)
data_raw <- temp

v <- NULL
w <- 0.1 # weight of the variance of the genomic region (1-w is the weight for the locus) 

# Calculates each locus' individual variance
for (i in chr) {
  current <- data_raw[which(data_raw[,1]==i),]
  vtemp <- as.data.frame(matrix(NA, nrow = nrow(current), ncol = 2))
  colnames(vtemp) <- c("v1", "v2")
  currentT <- t(current[,c(3:10)])
  
  for (j in c(1:nrow(current))) {
    vtemp[j,1] <- var(currentT[c(1:4),j])
    vtemp[j,2] <- var(currentT[c(5:8),j])
  }
  
  v <- rbind(v, vtemp)
}

# Calculate means for C and HR for each locus
data_rawc <- data_raw[,c(3:6)] 
data_rawHR <- data_raw[,c(7:10)]

data_rawc$x1 <- rowMeans(data_rawc)
data_rawHR$x2 <- rowMeans(data_rawHR)

accum <- cbind(data_raw[,c(1,2)], data_rawc[,5], data_rawHR[,5], v)

# data contains: "Chr", "POS", "x1" (C AF means), "x2" (HR AF means), "v1" (C AF variance), and "v2" (HR AF variance)
data <- accum[order(accum[,1],accum[,2]),]

data$vv <- NA # Will be the variance means for the region
final <- NULL

# Calculates sliding window variance
for (i in chr) {
  current <- data[which(data[,1]==i),]
  
  for (j in 1:nrow(current)) {
    subcurrent <- current[which(current[,2]>=current[j,2]-1000000 & current[,2]<=current[j,2]+1000000),]
    substack <- stack(subcurrent[,c(5,6)])
    current[j,7] <- mean(substack[,1])
  }
  final <- rbind(final, current)
}

results <- NULL

# The sliding window variance is written out as an intermediate file 
setwd(dir_out)
write.csv(final, paste("WRF_gen", generation, "_sliding_window_variance.csv", sep = ""), row.names = FALSE)

# This calculates the WRF using the above window variance
for (i in chr) {
  current <- final[which(final[,1]==i),]
  results1 <- as.data.frame(matrix(NA, nrow = nrow(current), ncol = 3))
  colnames(results1) <- c("regF", "regF_P", "DF") 
  
  for (j in c(1:nrow(current))) {
    results1[j,1] <- ((current[j,3] - current[j,4])^2)/(((1-w) * (current[j,5]/4 + current[j,6]/4)) + (2 * w * current[j,7])/4) 
    results1[j,3] <- (((current[j,5]/4)+(current[j,6]/4))**2)/((current[j,5]**2)/48+(current[j,6]**2)/48) # only need for different sample size
    results1[j,2] <- pf(results1[j,1], 1, results1[j,3], lower.tail = FALSE)
  }
  
  results <- rbind(results, results1)
  print(i)
}

final <- cbind(final, results)
colnames(final)[3] <- "C_means"
colnames(final)[4] <- "HR_means"

write.csv(final, paste("WRF_gen", generation, "_results.csv", sep = ""), row.names = FALSE)
