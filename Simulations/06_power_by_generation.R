# This code was written for: 
# David A Hillis, Liran Yadgary, George M. Weinstock, Fernando Pardo-Manuel
# de Villena, Daniel Pomp, and Theodore Garland Jr. “Large Changes in Detected Selection Signatures After a
# Selection Limit in Mice Bred for Voluntary Wheel-Running Behavior.” PLOS ONE, 2024.

# This calculates power by generation for the results of 04_calculate_power_for_multiple_iterations.R 

results_dir <- "/home/davidhillis/Simulation_trial/Results/"
setwd(results_dir)
letters <- LETTERS[c(1)] 

data <- NULL
for (i in letters) {data <- rbind(data, read.csv(paste("2096_Loci_multi_power", i, ".csv", sep = ""), stringsAsFactors = FALSE))}
datacon <- data[which(data[,2]=="constraint"),]
datauncon <- data[which(data[,2]=="unconstrained"),]

iter <- c(seq(0,20,5), 22, seq(25,60,5), 61)

final <- as.data.frame(matrix(NA, ncol = 5, nrow = 15))
final <- cbind(as.data.frame(iter), final)
colnames(final) <- c("generation", "Unconstrained_power", "p_difnextgen_uncon", "Constrained_power", "p_difnextgen_con", "p_Con_vs_uncon")

# Statistical tests for a difference in powere between generations
for (i in 1:nrow(final)) {
  final[i,2] <- mean(datauncon[,(i+2)])
  if (i!=nrow(final)) {test <- t.test(datauncon[,(i+2)], datauncon[,(i+3)])
  final[i,3] <- test$p.value}
  
  final[i,4] <- mean(datacon[,(i+2)])
  if (i!=nrow(final)) {test <- t.test(datacon[,(i+2)], datacon[,(i+3)])
  final[i,5] <- test$p.value}
  
  test <- t.test(datauncon[,(i+2)], datacon[,(i+2)])
  final[i,6] <- test$p.value
}

write.csv(final, "Power_comparisons_by_generation_and_constraint.csv", row.names = FALSE)

# Scatterplot of Generation 20
# The generation can be changed by changing the column referenced in datauncon and datacon
uncon <- cbind(as.data.frame(rep(1,100)), datauncon[,7])
colnames(uncon) <- c("model", "power")
con <- cbind(as.data.frame(rep(2,100)), datacon[,7])
colnames(con) <- c("model", "power")
newdata <- rbind(uncon, con)
plot(newdata[,1], newdata[,2])
t.test(uncon[,2], con[,2])

unconmean <- mean(uncon[,2])
conmean <- mean(con[,2])

write.csv(newdata, "gen20_power_comparisons_unconstrained1_constrained2.csv", row.names = FALSE)
jpeg(filename = paste("generation20_power_scatterplot.jpeg", sep = ""), width = 600, height = 3000, quality = 500)
par(mar = c(5, 6, 4, 2) + 0.1)
plot(jitter(newdata[,1]), newdata[,2], lwd = 5, axes = FALSE, cex.main = 5, cex.lab = 4, cex = 4, main = "", xlab = "Model", ylab = "Power")
lines(c(1,2), c(unconmean, conmean), lwd = 4, col = "red")
dev.off()


