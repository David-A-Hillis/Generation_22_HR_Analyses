# This code was written for: 
# David A Hillis, Liran Yadgary, George M. Weinstock, Fernando Pardo-Manuel
# de Villena, Daniel Pomp, and Theodore Garland Jr. “Large Changes in Detected Selection Signatures After a
# Selection Limit in Mice Bred for Voluntary Wheel-Running Behavior.” PLOS ONE, 2024.

# This performs simulations of allelic response to selection given various realistic factors

data_dir <- "/home/davidhillis/Simulation_trial/Data/" # identify working directory to place data
total_replications <- 10 # Total number of repeats for the simulations to conduct (each can take 6+ hours)
unconstrained_replications <- 5 # How many of the total replications to be unconstrained
N <- 20 # number of mice parenting each generation (half males and half females)
# Each of the 10 pairs will produce 1 male and 1 female for controls and "sel" males and "sel" females for HR as shown in next line
sel <- 5 # Number of mice being selected from (for each sex)
c_sel <- 2 # Number of CONTROL mice produced for each sex (should be 1 unless doing diagnostics)
loci <- 2096 # number of loci influencing the trait (total wheel running) NOTE Will create problems with kurtosis if exceeds 4000
start_af <- 0.5 # starting allele frequency for loci (0.25 was picked for maintaining running at about 6000 for controls and allows for a total of 21000)
gen <- 61 # Number of generations
constraint1 <- 10000 # number of revolutions that cannot be exceeded 
Ve_constraint <- 1750 # Environmental multiplier for constraint
max_running <- 50000 # Highest running allowed for a mouse, running levels above this are replaced with this value
Min_running <- 100 # Lowest running allowed for a mouse, running levels lower this are replaced with this value
gene_mod <- 1.3 # Multiplier for the genetics variance (larger numbers increases running variance)
Mean_start <- 4570 # Starting population mean (1527.4 according to my calculations based on average gain from dominance [complete, 1/4], 4400 if no dominance)
Ve <- 2100 # Set the multiplier for the environmental variance, normally 500
dominance <- 0 # amount of dominance 0 to 2 (e.g., 1 for dominant heterozygotes to have 50% genetic effect as homozygous dominant)
dom_proportion <- 1/4 # proportion of loci to be affected by dominance, should be listed as 1/x (e.g., 1/5 for 20% dominance rate)
season <- c(0.769, 1, 1.3, 1) # Vector representing running multiplier in each season c(Summer, Fall, Winter, Spring), not used for gen 0
s_loci <- 100 # Number of loci to influence seasonal response
s_start_af <- 0.5 # Starting allele frequency for seasonal alleles 
s_epistatis <- 0 # Set to 1 in order to include selection on seasonal response
seasonalVar <- 0 # set equal to 1 if you want seasonal variation by locus
con_loci <- 100 # Number of loci to influence seasonal response
con_start_af <- 0.5 # Starting allele frequency for seasonal alleles 
con_epistatis <- 1 # Set to 1 in order to include selection on seasonal response
write_breed <- 1 # set to 1 if you want to write out the breed order or 0 if not (needed for heritability)
# (i.e., maximum of 21000 occurs when all alleles = 1, and 1000 is the base value; note that base value sets the min)

# Seasonal variation by locus (not used if seasonalVar doesn't equal 1)
# summer <- rep(c(1.7,0.3,1,1, 1,1.7,0.3,1, 1,1,1.7,0.3, 0.3,1,1,1.7),131)
# winter <- rep(c(0.3,1.7,1,1, 1,0.3,1.7,1, 1,1,0.3,1.7, 1.7,1,1,0.3),131)

# Create Weight of Constraint Loci
if (con_epistatis==1) {conweight <- c(rep(1,48), rep(4,24), rep(11,12), rep(36,8), rep(101,5), rep(306,3))}
save_con_0 <- NULL
save_con_C <- NULL
save_con_HR <- NULL

set.seed(42) # seed 43 for B, and seed 44 for C

for (z in 1:total_replications) { # Change for the number of iterations to be completed
  setwd(data_dir)
  if (z<=unconstrained_replications) {constraint <- max_running}
  else {constraint <- constraint1}
  
  running <- as.data.frame(matrix(NA, nrow = 0, ncol = (N+1)))
  colrun <- c("line")
  for (i in 1:N) {colrun <- c(colrun, paste("m", i, sep = ""))}
  colnames(running) <- colrun
  
  # Implement kurtosis (allelic weights)
  gweight <- c(rep(0.4, 720), rep(0.8, 480), rep(1.6, 312), rep(3.2, 216), rep(6.4, 144),
               rep(12.8, 96), rep(25.6, 60), rep(51.2, 36), rep(102.4, 24), rep(204.8, 8))
  
  # Create Weight of seasonal loci (Total to 0.768)
  if (s_epistatis==1) {sweight <- c(rep(0.001,48), rep(0.003,24), rep(0.005,12), rep(0.014,8), rep(0.04,5), rep(0.120,3))}
  
  # control mice
  wholeRunning <- as.data.frame(matrix(NA, nrow = 0, ncol = N*c_sel+1))
  cnam <- c("line")
  for (j in 1:(N/2)) {for (k in 1:2) {for (m in 1:c_sel) {cnam <- c(cnam, paste("Fam", j, "_Sex", k, "_m", m, sep = ""))}}}
  colnames(wholeRunning) <- cnam
  
  for (i in 1:4) {
    # Create starting matrix for each line
    if (con_epistatis==1) {
      con_alleles <- as.data.frame(matrix(rbinom((2*N*con_loci),1,con_start_af), ncol = (2*N), nrow = con_loci))
      for (j in 1:nrow(con_alleles)) {for (k in 1:ncol(con_alleles)) {if (con_alleles[j,k]==0) {con_alleles[j,k] <- -1}}}
      for (j in 1:nrow(con_alleles)) {for (k in 1:ncol(con_alleles)) {con_alleles[j,k] <- con_alleles[j,k]*conweight[j]}}
    }
    
    if (s_epistatis==1) {
      s_alleles <- as.data.frame(matrix(rbinom((2*N*s_loci),1,s_start_af), ncol = (2*N), nrow = s_loci)) # create starting alleles
      for (j in 1:nrow(s_alleles)) {for (k in 1:ncol(s_alleles)) {if (s_alleles[j,k]==0) {s_alleles[j,k] <- -1}}}
      for (j in 1:nrow(s_alleles)) {for (k in 1:ncol(s_alleles)) {s_alleles[j,k] <- s_alleles[j,k]*sweight[j]}}
    }
    
    alleles <- as.data.frame(matrix(rbinom((2*N*loci),1,start_af), ncol = (2*N), nrow = loci)) # create starting alleles
    for (j in 1:nrow(alleles)) {for (k in 1:ncol(alleles)) {if (alleles[j,k]==0) {alleles[j,k] <- -1}}}
    for (j in 1:nrow(alleles)) {for (k in 1:ncol(alleles)) {alleles[j,k] <- alleles[j,k]*gweight[j]}}
    
    ### Calculate and save running 
    if (dominance == 0) {af <- colSums(alleles)}
    else {
      af <- c()
      for (j in 1:N) {
        af_temp <- alleles[,((((j-1)*2)+1):(((j-1)*2)+2))]
        af_temp$sums <- rowSums(af_temp)
        for (k in 1:loci) {if (k%%(1/dom_proportion)==0 & af_temp[k,3]==0) {af_temp[k,3] <- dominance*max(af_temp[k,])}}
        af <- c(af, sum(af_temp[,3]), 0)
      }
    }
    
    distance <- NULL
    # Add seasonal epistasis to base population running levels
    if (s_epistatis==1) {
      s_af <- colSums(s_alleles)
      for (j in 1:N) {distance <- c(distance,
                                    Mean_start+(rnorm(1)*Ve)+(s_af[((j-1)*2)+1]+s_af[((j-1)*2)+2]+1)*(gene_mod)*(af[((j-1)*2)+1]+af[((j-1)*2)+2]))}
      write.csv(s_alleles, paste("season_alleles_line", i, "_", s_loci, "loci_gen0_runA", z, ".csv", sep = ""), row.names = FALSE)
    }
    else {for (j in 1:N) {distance <- c(distance, Mean_start+(rnorm(1)*Ve)+(gene_mod)*(af[((j-1)*2)+1]+af[((j-1)*2)+2]))}} # No seasonal variation for gen 0
    distance[which(distance<Min_running)] <- Min_running
    
    # Change set max distance for base generation
    if (con_epistatis==1) {
      #      if (s_epistatis==1) {con_af <- colSums(con_alleles)+(constraint*season[(k%%4)+1])} # Removed to untangle con and season
      con_af <- colSums(con_alleles)
      con_temp <- c()
      for (j in 1:length(distance)) {
        temp_vec <- (rnorm(1)*Ve_constraint)
        con_temp <- c(con_temp, con_af[((j-1)*2)+1]+con_af[((j-1)*2)+2]+constraint+temp_vec)
        if (distance[j]>(con_af[((j-1)*2)+1]+con_af[((j-1)*2)+2]+constraint+temp_vec)) 
        {distance[j] <- (con_af[((j-1)*2)+1]+con_af[((j-1)*2)+2]+constraint+temp_vec)}}
      save_con_0 <- rbind(save_con_0, t(as.data.frame(c(i, 0, con_temp))))
    } 
    else {distance[which(distance>constraint)] <- constraint}
    
    distance <- t(as.data.frame(c(i, distance)))
    colnames(distance) <- colnames(running)
    running <- rbind(running, distance)
    write.csv(alleles, paste("constraint_simulation_line", i, "_", loci, "loci_gen0_runA", z, ".csv", sep = ""), row.names = FALSE)
    
    #start drift simulation
    for (k in 1:(gen)) {
      
      # determine pairs
      b1 <- c(1:(N/2))
      this <- 1
      while (this==1) {
        temp_parent <- sample(b1,N/2) # Realistically, the larger N gets, the more likely this will have to run a while to complete.
        breed_order <- cbind(b1, as.data.frame(temp_parent))
        if (sum(breed_order[,1]==breed_order[,2])>0) {}
        else {this <- 0}
      }
      if (write_breed==1) {write.csv(breed_order, paste("breedOrder_line", i, "_", loci, "loci_gen", k,"_runA", z, ".csv", sep = ""), row.names = FALSE)}
      
      # randomly pass alleles
      allele_temp <- as.data.frame(matrix(NA, ncol = (2*N*c_sel), nrow = loci))
      
      for (j in 1:(N/2)) { # each breeding pair
        for (m in 1:loci) { # each locus
          for (n in 1:(c_sel*2)) { # each offspring
            allele_temp[m,(((j-1)*(c_sel*4))+(((n-1)*2)+1))] <- as.numeric(sample(alleles[m,c(((breed_order[j,1]-1)*4)+1, ((breed_order[j,1]-1)*4)+2)], 1))
            allele_temp[m,(((j-1)*(c_sel*4))+(((n-1)*2)+2))] <- as.numeric(sample(alleles[m,c(((breed_order[j,2]-1)*4)+3, ((breed_order[j,2]-1)*4)+4)], 1))
          }
        }
      }      
      alleles <- as.data.frame(matrix(NA, nrow = loci, ncol = 0))
      for (j in 1:N) {alleles <- cbind(alleles, allele_temp[,c((((j-1)*(2*c_sel))+1), (((j-1)*(2*c_sel))+2))])}
      
      # Randomly pass Constraint alleles
      if (con_epistatis==1) {
        con_allele_temp <- as.data.frame(matrix(NA, ncol = (2*N*c_sel), nrow = con_loci))
        
        for (j in 1:(N/2)) { # each breeding pair
          for (m in 1:con_loci) { # each locus
            for (n in 1:(c_sel*2)) { # each offspring
              con_allele_temp[m,(((j-1)*(c_sel*4))+(((n-1)*2)+1))] <- as.numeric(sample(con_alleles[m,c(((breed_order[j,1]-1)*4)+1, ((breed_order[j,1]-1)*4)+2)], 1))
              con_allele_temp[m,(((j-1)*(c_sel*4))+(((n-1)*2)+2))] <- as.numeric(sample(con_alleles[m,c(((breed_order[j,2]-1)*4)+3, ((breed_order[j,2]-1)*4)+4)], 1))
            }
          }
        }      
        con_alleles <- as.data.frame(matrix(NA, nrow = con_loci, ncol = 0))
        for (j in 1:N) {con_alleles <- cbind(con_alleles, con_allele_temp[,c((((j-1)*(2*c_sel))+1), (((j-1)*(2*c_sel))+2))])}
      }
      
      # Randomly pass seasonal alleles
      if (s_epistatis==1) {
        s_allele_temp <- as.data.frame(matrix(NA, ncol = (2*N*c_sel), nrow = s_loci))
        
        for (j in 1:(N/2)) { # each breeding pair
          for (m in 1:s_loci) { # each locus
            for (n in 1:(c_sel*2)) { # each offspring
              s_allele_temp[m,(((j-1)*(c_sel*4))+(((n-1)*2)+1))] <- as.numeric(sample(s_alleles[m,c(((breed_order[j,1]-1)*4)+1, ((breed_order[j,1]-1)*4)+2)], 1))
              s_allele_temp[m,(((j-1)*(c_sel*4))+(((n-1)*2)+2))] <- as.numeric(sample(s_alleles[m,c(((breed_order[j,2]-1)*4)+3, ((breed_order[j,2]-1)*4)+4)], 1))
            }
          }
        }      
        s_alleles <- as.data.frame(matrix(NA, nrow = s_loci, ncol = 0))
        for (j in 1:N) {s_alleles <- cbind(s_alleles, s_allele_temp[,c((((j-1)*(2*c_sel))+1), (((j-1)*(2*c_sel))+2))])}
      }
      
      # Add individual locus seasonal variation
      #      if (seasonalVar==1 & k%%4==2) {for (j in (1:nrow(allele_temp))) {allele_temp[j,] <- allele_temp[j,]*winter[j]}} # Have not been used
      #      else if (seasonalVar==1 & k%%4==0) {for (j in (1:nrow(allele_temp))) {allele_temp[j,] <- allele_temp[j,]*summer[j]}} # Have not been used
      
      ### Calculate and save running 
      if (dominance == 0) {af_temp1 <- colSums(allele_temp)}
      else {
        af_temp1 <- c()
        for (j in 1:(N*c_sel)) {
          af_temp <- allele_temp[,((((j-1)*2)+1):(((j-1)*2)+2))]
          af_temp$sums <- rowSums(af_temp)
          for (l in 1:loci) {if (l%%(1/dom_proportion)==0 & af_temp[l,3]==0) {af_temp[l,3] <- dominance*max(af_temp[l,])}}
          af_temp1 <- c(af_temp1, sum(af_temp[,3]), 0)
        }
      }
      
      distance <- NULL
      
      # Add seasonal affect to next generation of offspring
      if (s_epistatis==1) {
        s_af_temp1 <- colSums(s_allele_temp)
        for (j in 1:(N*c_sel)) {distance <- c(distance, season[(k%%4)+1]*(Mean_start+(rnorm(1)*Ve)+(s_af_temp1[((j-1)*2)+1]+s_af_temp1[((j-1)*2)+2]+1)*(gene_mod)*(af_temp1[((j-1)*2)+1]+af_temp1[((j-1)*2)+2])))}
        if (k%%5==0 | k==22 | k==61) {write.csv(s_alleles, paste("season_alleles_line", i, "_", loci, "loci_gen", k, "_runA", z, ".csv", sep = ""), row.names = FALSE)}
      }
      else {for (j in 1:(N*c_sel)) {distance <- c(distance, season[(k%%4)+1]*(Mean_start+(rnorm(1)*Ve)+(gene_mod)*(af_temp1[((j-1)*2)+1]+af_temp1[((j-1)*2)+2])))}}
      distance[which(distance<Min_running)] <- Min_running
      
      # Change set max distance
      if (con_epistatis==1) {
        #        if (s_epistatis==1) {con_af_temp1 <- colSums(con_allele_temp)+(constraint*season[(k%%4)+1])} # Decouple season and multiplier genes
        con_af_temp1 <- colSums(con_allele_temp)
        con_temp <- c()
        for (j in 1:length(distance)) {
          temp_vec <- (rnorm(1)*Ve_constraint)
          con_temp <- c(con_temp, season[(k%%4)+1]*con_af_temp1[((j-1)*2)+1]+con_af_temp1[((j-1)*2)+2]+constraint+temp_vec)
          if (distance[j]>season[(k%%4)+1]*(con_af_temp1[((j-1)*2)+1]+con_af_temp1[((j-1)*2)+2]+constraint+temp_vec)) 
          {distance[j] <- season[(k%%4)+1]*(con_af_temp1[((j-1)*2)+1]+con_af_temp1[((j-1)*2)+2]+constraint+temp_vec)}}
        if (k%%5==0 | k==22 | k==61) {write.csv(con_allele_temp, paste("Constraint_alleles_line", i, "_", loci, "loci_gen", k, "_runA", z, ".csv", sep = ""), row.names = FALSE)}
        save_con_C <- rbind(save_con_C, t(as.data.frame(c(i, k, con_temp[seq(1,(N*c_sel-1),c_sel)]))))
      } 
      else {distance[which(distance>constraint)] <- constraint}
      
      distance <- t(as.data.frame(c(i, distance)))
      colnames(distance) <- colnames(wholeRunning)
      wholeRunning <- rbind(wholeRunning, distance)
      
      distance <- t(as.data.frame(distance[,c(1,seq(2,(2+((N-1)*c_sel)),c_sel))]))
      colnames(distance) <- colnames(running)
      running <- rbind(running, distance)
      
      if (k%%5==0 | k==22 | k==61) {write.csv(alleles, 
                                              paste("constraint_simulation_line", i, "_", loci, "loci_gen", k, "_runA", z, ".csv", sep = ""), row.names = FALSE)}
    }
    write.csv(alleles, paste("constraint_simulation_line", i, "_", loci, "loci_complete_C", "_runA", z, ".csv", sep = ""), row.names = FALSE)
  }
  write.csv(wholeRunning, paste("wholeRunning_C_", loci, "loci_50AF_v77_constraint", constraint, "_runA", z, ".csv", sep = ""), row.names = FALSE)
  
  ########################
  ### High-Runner Mice ###
  ########################
  
  wholeRunning <- as.data.frame(matrix(NA, nrow = 0, ncol = N*sel+1))
  cnam <- c("line")
  for (j in 1:(N/2)) {for (k in 1:2) {for (m in 1:sel) {cnam <- c(cnam, paste("Fam", j, "_Sex", k, "_m", m, sep = ""))}}}
  colnames(wholeRunning) <- cnam
  
  for (i in 5:8) {
    # Create starting matrix for each line
    if (con_epistatis==1) {
      con_alleles <- as.data.frame(matrix(rbinom((2*N*con_loci),1,con_start_af), ncol = (2*N), nrow = con_loci))
      for (j in 1:nrow(con_alleles)) {for (k in 1:ncol(con_alleles)) {if (con_alleles[j,k]==0) {con_alleles[j,k] <- -1}}}
      for (j in 1:nrow(con_alleles)) {for (k in 1:ncol(con_alleles)) {con_alleles[j,k] <- con_alleles[j,k]*conweight[j]}}
    }
    
    if (s_epistatis==1) {
      s_alleles <- as.data.frame(matrix(rbinom((2*N*s_loci),1,s_start_af), ncol = (2*N), nrow = s_loci)) # create starting alleles
      for (j in 1:nrow(s_alleles)) {for (k in 1:ncol(s_alleles)) {if (s_alleles[j,k]==0) {s_alleles[j,k] <- -1}}}
      for (j in 1:nrow(s_alleles)) {for (k in 1:ncol(s_alleles)) {s_alleles[j,k] <- s_alleles[j,k]*sweight[j]}}
    }
    
    alleles <- as.data.frame(matrix(rbinom((2*N*loci),1,start_af), ncol = (2*N), nrow = loci)) # create starting alleles
    for (j in 1:nrow(alleles)) {for (k in 1:ncol(alleles)) {if (alleles[j,k]==0) {alleles[j,k] <- -1}}}
    for (j in 1:nrow(alleles)) {for (k in 1:ncol(alleles)) {alleles[j,k] <- alleles[j,k]*gweight[j]}}
    
    ### Calculate and save running
    if (dominance == 0) {af <- colSums(alleles)}
    else {
      af <- c()
      for (j in 1:N) {
        af_temp <- alleles[,((((j-1)*2)+1):(((j-1)*2)+2))]
        af_temp$sums <- rowSums(af_temp)
        for (l in 1:loci) {if (l%%(1/dom_proportion)==0 & af_temp[l,3]==0) {af_temp[l,3] <- dominance*max(af_temp[l,])}}
        af <- c(af, sum(af_temp[,3]), 0)
      }
    }
    
    distance <- NULL
    
    # Add seasonal epistasis to base population running levels
    if (s_epistatis==1) {
      s_af <- colSums(s_alleles)
      for (j in 1:N) {distance <- c(distance,
                                    Mean_start+(rnorm(1)*Ve)+(s_af[((j-1)*2)+1]+s_af[((j-1)*2)+2]+1)*(gene_mod)*(af[((j-1)*2)+1]+af[((j-1)*2)+2]))}
      write.csv(s_alleles, paste("season_alleles_line", i, "_", s_loci, "loci_gen0_runA", z, ".csv", sep = ""), row.names = FALSE)
    }
    else {for (j in 1:N) {distance <- c(distance, Mean_start+(rnorm(1)*Ve)+(gene_mod)*(af[((j-1)*2)+1]+af[((j-1)*2)+2]))}}
    distance[which(distance<Min_running)] <- Min_running
    
    # Change set max distance for base generation
    if (con_epistatis==1) {
      #      if (s_epistatis==1) {con_af <- colSums(con_alleles)+(constraint*season[(k%%4)+1])} # Detangle
      con_af <- colSums(con_alleles)
      con_temp <- c()
      for (j in 1:length(distance)) {
        temp_vec <- (rnorm(1)*Ve_constraint)
        con_temp <- c(con_temp, con_af[((j-1)*2)+1]+con_af[((j-1)*2)+2]+constraint+temp_vec)
        if (distance[j]>(con_af[((j-1)*2)+1]+con_af[((j-1)*2)+2]+constraint+temp_vec)) 
        {distance[j] <- (con_af[((j-1)*2)+1]+con_af[((j-1)*2)+2]+constraint+temp_vec)}}
      save_con_0 <- rbind(save_con_0, t(as.data.frame(c(i, 0, con_temp))))
    }
    else {distance[which(distance>constraint)] <- constraint}
    
    distance <- t(as.data.frame(c(i, distance)))
    colnames(distance) <- colnames(running)
    running <- rbind(running, distance)
    write.csv(alleles, paste("constraint_simulation_line", i, "_", loci, "loci_gen0_runA", z, ".csv", sep = ""), row.names = FALSE)
    
    # File for saving out all offspring running
    wholeRunning_temp <- as.data.frame(matrix(NA, nrow = 0, ncol = N*sel+2))
    colnames(wholeRunning_temp) <- cnam
    ###
    
    #start drift simulation
    for (k in 1:gen) {
      
      # determine pairs
      b1 <- c(1:(N/2))
      this <- 1
      while (this==1) {
        temp_parent <- sample(b1,N/2)
        breed_order <- cbind(b1, as.data.frame(temp_parent))
        if (sum(breed_order[,1]==breed_order[,2])>0) {}
        else {this <- 0}
      }
      if (write_breed==1) {write.csv(breed_order, paste("breedOrder_line", i, "_", loci, "loci_gen", k,"_runA", z, ".csv", sep = ""), row.names = FALSE)}
      
      # randomly pass alleles
      allele_temp <- as.data.frame(matrix(NA, ncol = (2*N*sel), nrow = loci))
      
      for (j in 1:(N/2)) { # each breeding pair
        for (m in 1:loci) { # each locus
          for (n in 1:(sel*2)) { # each offspring
            allele_temp[m,(((j-1)*(sel*4))+(((n-1)*2)+1))] <- as.numeric(sample(alleles[m,c(((breed_order[j,1]-1)*4)+1, ((breed_order[j,1]-1)*4)+2)], 1))
            allele_temp[m,(((j-1)*(sel*4))+(((n-1)*2)+2))] <- as.numeric(sample(alleles[m,c(((breed_order[j,2]-1)*4)+3, ((breed_order[j,2]-1)*4)+4)], 1))
          }
        }
      }
      
      # Randomly pass Constraint alleles
      if (con_epistatis==1) {
        con_allele_temp <- as.data.frame(matrix(NA, ncol = (2*N*sel), nrow = con_loci))
        
        for (j in 1:(N/2)) { # each breeding pair
          for (m in 1:con_loci) { # each locus
            for (n in 1:(sel*2)) { # each offspring
              con_allele_temp[m,(((j-1)*(sel*4))+(((n-1)*2)+1))] <- as.numeric(sample(con_alleles[m,c(((breed_order[j,1]-1)*4)+1, ((breed_order[j,1]-1)*4)+2)], 1))
              con_allele_temp[m,(((j-1)*(sel*4))+(((n-1)*2)+2))] <- as.numeric(sample(con_alleles[m,c(((breed_order[j,2]-1)*4)+3, ((breed_order[j,2]-1)*4)+4)], 1))
            }
          }
        }
      }
      
      # Randomly pass seasonal alleles
      if (s_epistatis==1) {
        s_allele_temp <- as.data.frame(matrix(NA, ncol = (2*N*sel), nrow = s_loci))
        
        for (j in 1:(N/2)) { # each breeding pair
          for (m in 1:s_loci) { # each locus
            for (n in 1:(sel*2)) { # each offspring
              s_allele_temp[m,(((j-1)*(sel*4))+(((n-1)*2)+1))] <- as.numeric(sample(s_alleles[m,c(((breed_order[j,1]-1)*4)+1, ((breed_order[j,1]-1)*4)+2)], 1))
              s_allele_temp[m,(((j-1)*(sel*4))+(((n-1)*2)+2))] <- as.numeric(sample(s_alleles[m,c(((breed_order[j,2]-1)*4)+3, ((breed_order[j,2]-1)*4)+4)], 1))
            }
          }
        }
      }
      
      # Add individual locus seasonal variation
      allele_season <- allele_temp
      #      if (seasonalVar==1 & k%%4==2) {for (j in (1:nrow(allele_season))) {allele_season[j,] <- allele_season[j,]*winter[j]}}
      #      else if (seasonalVar==1 & k%%4==0) {for (j in (1:nrow(allele_season))) {allele_season[j,] <- allele_season[j,]*summer[j]}}
      
      ### Calculate running for everyone
      if (dominance == 0) {af_temp1 <- colSums(allele_season)}
      else {
        af_temp1 <- c()
        for (j in 1:(N*sel)) {
          af_temp <- allele_season[,((((j-1)*2)+1):(((j-1)*2)+2))]
          af_temp$sums <- rowSums(af_temp)
          for (l in 1:loci) {if (l%%(1/dom_proportion)==0 & af_temp[l,3]==0) {af_temp[l,3] <- dominance*max(af_temp[l,])}}
          af_temp1 <- c(af_temp1, sum(af_temp[,3]), 0)
        }
      }
      
      distance_temp <- c()
      
      # Add seasonal affect to next generation of offspring
      if (s_epistatis==1) {
        s_af_temp1 <- colSums(s_allele_temp)
        for (j in 1:(N*sel)) {distance_temp <- c(distance_temp, season[(k%%4)+1]*(Mean_start+(rnorm(1)*Ve)+(s_af_temp1[((j-1)*2)+1]+s_af_temp1[((j-1)*2)+2]+1)*(gene_mod)*(af_temp1[((j-1)*2)+1]+af_temp1[((j-1)*2)+2])))}
        if (k%%5==0 | k==22 | k==61) {write.csv(s_alleles, paste("season_alleles_line", i, "_", loci, "loci_gen", k, "_runA", z, ".csv", sep = ""), row.names = FALSE)}
      }
      else {for (j in 1:(N*sel)) {distance_temp <- c(distance_temp, season[(k%%4)+1]*(Mean_start+(rnorm(1)*Ve)+(gene_mod)*(af_temp1[((j-1)*2)+1]+af_temp1[((j-1)*2)+2])))}}
      distance_temp[which(distance_temp<Min_running)] <- Min_running
      
      # Change set max distance
      if (con_epistatis==1) {
        #        if (s_epistatis==1) {con_af_temp1 <- colSums(con_allele_temp)+(constraint*season[(k%%4)+1])} # Decouple
        con_af_temp1 <- colSums(con_allele_temp)
        con_temp <- c()
        for (j in 1:length(distance_temp)) {
          temp_vec <- (rnorm(1)*Ve_constraint)
          con_temp <- c(con_temp, season[(k%%4)+1]*con_af_temp1[((j-1)*2)+1]+con_af_temp1[((j-1)*2)+2]+constraint+temp_vec)
          if (distance_temp[j]>season[(k%%4)+1]*(con_af_temp1[((j-1)*2)+1]+con_af_temp1[((j-1)*2)+2]+constraint+temp_vec)) 
          {distance_temp[j] <- season[(k%%4)+1]*(con_af_temp1[((j-1)*2)+1]+con_af_temp1[((j-1)*2)+2]+constraint+temp_vec)}}
      }
      else {distance_temp[which(distance_temp>constraint)] <- constraint}
      
      # select breeders
      af2 <- distance_temp
      keep <- c()
      for (j in 1:N) {
        if (length(which(af2[(((j-1)*sel)+1):(((j-1)*sel)+sel)]<=Min_running))==sel) {
          keep <- c(keep, (as.numeric(sample(which(af2[(((j-1)*sel)+1):(((j-1)*sel)+sel)]<=Min_running), 1))+((j-1)*sel)))}
        else if (con_epistatis==1) {
          temp_max <- which(af2[(((j-1)*sel)+1):(((j-1)*sel)+sel)]==max(af2[(((j-1)*sel)+1):(((j-1)*sel)+sel)]))
          keep <- c(keep, (temp_max[1]+((j-1)*sel)))}
        else if (length(which(af2[(((j-1)*sel)+1):(((j-1)*sel)+sel)]>=constraint))<2) {
          temp_max <- which(af2[(((j-1)*sel)+1):(((j-1)*sel)+sel)]==max(af2[(((j-1)*sel)+1):(((j-1)*sel)+sel)]))
          keep <- c(keep, (temp_max+((j-1)*sel)))}
        else {keep <- c(keep, (as.numeric(sample(which(af2[(((j-1)*sel)+1):(((j-1)*sel)+sel)]>=constraint), 1))+((j-1)*sel)))}
      }
      
      alleles <- as.data.frame(matrix(NA, ncol = 0, nrow = loci))
      for (j in keep) {alleles <- cbind(alleles, allele_temp[,c((((j-1)*2)+1), (((j-1)*2)+2))])}
      
      # Apply selection to Constraint alleles
      if (con_epistatis==1) {
        con_alleles <- as.data.frame(matrix(NA, ncol = 0, nrow = con_loci))
        for (j in keep) {con_alleles <- cbind(con_alleles, con_allele_temp[,c((((j-1)*2)+1), (((j-1)*2)+2))])}
        con_af <- con_temp[keep]
        if (k%%5==0 | k==22 | k==61) {write.csv(con_allele_temp, paste("Constraint_alleles_line", i, "_", loci, "loci_gen", k, "_runA", z, ".csv", sep = ""), row.names = FALSE)}
        save_con_HR <- rbind(save_con_HR, t(as.data.frame(c(i, k, con_af))))
      }
      
      # Apply selection to seasonal alleles
      if (s_epistatis==1) {
        s_alleles <- as.data.frame(matrix(NA, ncol = 0, nrow = s_loci))
        for (j in keep) {s_alleles <- cbind(s_alleles, s_allele_temp[,c((((j-1)*2)+1), (((j-1)*2)+2))])}
        s_af <- colSums(s_alleles)
        if (k%%5==0 | k==22 | k==61) {write.csv(s_allele_temp, paste("season_alleles_line", i, "_", loci, "loci_gen", k, "_runA", z, ".csv", sep = ""), row.names = FALSE)}
      }
      
      # Save wheel running of breeders
      af <- colSums(alleles)
      distance <- distance_temp[keep]
      distance <- t(as.data.frame(c(i, distance)))
      colnames(distance) <- colnames(running)
      running <- rbind(running, distance)
      
      # Save all running
      distance_temp <- t(as.data.frame(c(i, distance_temp)))
      colnames(distance_temp) <- colnames(wholeRunning)
      wholeRunning_temp <- rbind(wholeRunning_temp, distance_temp)
      ###
      
      # First HR male and female for each family saved independent of running
      if (k%%5==0 | k==22 | k==61) {write.csv(allele_temp,
                                              paste("constraint_simulation_line", i, "_", loci, "loci_gen", k, "_runA", z, ".csv", sep = ""), row.names = FALSE)}
    }
    write.csv(alleles, paste("constraint_simulation_line", i, "_", loci, "loci_complete_HR", "_runA", z, ".csv", sep = ""), row.names = FALSE)
    wholeRunning <- rbind(wholeRunning, wholeRunning_temp)
  }
  
  write.csv(save_con_0, paste("constraint_levels_0_runA", z, ".csv", sep = ""), row.names = FALSE)
  write.csv(save_con_C, paste("constraint_levels_C_runA", z, ".csv", sep = ""), row.names = FALSE)
  write.csv(save_con_HR, paste("constraint_levels_HR_runA", z, ".csv", sep = ""), row.names = FALSE)
  write.csv(running, paste("constraint_simulation_", loci, "loci_50AF_v77_constraint", constraint, "_runA", z, ".csv", sep = ""), row.names = FALSE)
  write.csv(wholeRunning, paste("wholeRunning_HR_", loci, "loci_50AF_v77_constraint", constraint, "_runA", z, ".csv", sep = ""), row.names = FALSE)
}