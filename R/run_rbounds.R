library(rbounds)

setwd("D:/Nextcloud/Causal_inference_project/R")

# Read features_table
features_table <- read.csv('../results/normal_samples_only/features_table.txt',header=TRUE,row.names=1)
features_table$currently_smoking <- as.logical(features_table[,"currently_smoking"])
counts_per_million <- read.csv('../results/normal_samples_only/counts_per_million_diff_genes.txt',header=TRUE,row.names=1)

###########################################################################################
# Perform genetic optimal matching and sensitivity analysis with hlsens
ATE_genetic <- data.frame()

for (gene in row.names(counts_per_million)){
  print(gene)
  features_table$counts_per_million <- as.list(counts_per_million[gene,])
  
  Y <- features_table$counts_per_million
  Tr <- features_table$currently_smoking
  
  #Now With GenMatch
  X <- cbind(features_table$age, features_table$gender, features_table$bmi)
  BalanceMat <- cbind(features_table$age, features_table$gender, features_table$bmi)
  
  #Genetic Weights
  gen1 <- GenMatch(Tr=Tr, X=X, BalanceMat=BalanceMat, estimand='ATE', pop.size=50, replace=FALSE, print=0)
  
  #Match
  mgen1 <- Match(Y=Y, Tr=Tr, X=X, Weight.matrix=gen1, replace=FALSE)
  
  #Sensitivity Tests
  hlsens_table <- hlsens(mgen1, Gamma=2.3, GammaInc=.1, pr=.1)
  
  ATE_genetic[gene,"ate"] <-hlsens_table$bounds[1,2]
  ATE_genetic[gene,"lower_bound"] <- tail(hlsens_table$bounds,n=1)[2]
  ATE_genetic[gene,"upper_bound"] <- tail(hlsens_table$bounds,n=1)[3]
}

write.table(ATE_genetic, '../results/normal_samples_only/ate_genetic_sensitivity.txt', sep=',')

###########################################################################################
# Load 1:1 ps matching and perform sensitivity analysis with hlsens

matching <- read.csv('../results/normal_samples_only/ps_matching_with_replacement_11.txt')
ATE_ps <- data.frame()

for (gene in row.names(counts_per_million)){
  print(gene)
  features_table$counts_per_million <- as.list(counts_per_million[gene,])
  
  Y0 <- c()
  Y1 <- c()
  for (row in row.names(matching)){
    Y0 <- c(Y0, features_table[c(toString(matching[row,"control_indices"])),"counts_per_million"])
    Y1 <- c(Y1, features_table[c(toString(matching[row,"treatment_indices"])),"counts_per_million"])
  }
  
  #Sensitivity Tests
  hlsens_table <- hlsens(as.numeric(Y1), as.numeric(Y0), Gamma=2.3, GammaInc=.1, pr=.1)
  
  ATE_ps[gene,"ate"] <-hlsens_table$bounds[1,2]
  ATE_ps[gene,"lower_bound"] <- tail(hlsens_table$bounds,n=1)[2]
  ATE_ps[gene,"upper_bound"] <- tail(hlsens_table$bounds,n=1)[3]
}

write.table(ATE_ps, '../results/normal_samples_only/ate_ps_sensitivity.txt', sep=',')

###########################################################################################
# Load 1:1 cov matching and perform sensitivity analysis with hlsens

matching <- read.csv('../results/normal_samples_only/cov_matching_with_replacement_11.txt')
ATE_cov <- data.frame()

for (gene in row.names(counts_per_million)){
  print(gene)
  features_table$counts_per_million <- as.list(counts_per_million[gene,])
  
  Y0 <- c()
  Y1 <- c()
  for (row in row.names(matching)){
    Y0 <- c(Y0, features_table[c(toString(matching[row,"control_indices"])),"counts_per_million"])
    Y1 <- c(Y1, features_table[c(toString(matching[row,"treatment_indices"])),"counts_per_million"])
  }
  
  #Sensitivity Tests
  hlsens_table <- hlsens(as.numeric(Y1), as.numeric(Y0), Gamma=2.3, GammaInc=.1, pr=.1)
  
  ATE_cov[gene,"ate"] <-hlsens_table$bounds[1,2]
  ATE_cov[gene,"lower_bound"] <- tail(hlsens_table$bounds,n=1)[2]
  ATE_cov[gene,"upper_bound"] <- tail(hlsens_table$bounds,n=1)[3]
}

write.table(ATE_cov, '../results/normal_samples_only/ate_cov_sensitivity.txt', sep=',')