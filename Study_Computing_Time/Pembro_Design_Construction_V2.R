# ===================================================================
# Comparisons of three adaptive randomization methods for
# continuous covariates and two treatments.
#
# Construction of designs using different strategies and parameter settings.
# 
# ===================================================================
library(randomizeR)
library(e1071) # For "moments" function.

source("adaptive_randomization.R")

max.sim <- 1000

# Set the parameters for all algorithms.
n.zero <- 8 
bc.prob <- 0.8

# Set parameters for Bertsimas et al's approach.
rho.p <- 6
upper.gamma.bound <- 4

# Set parameters for Pocock and Simons.
n.categories <- 3

# Develop simulation protocol.
p <- 3 # Number of covariates.
my.data <- read.table("data/Pembrolizumab.txt", header = T)
my.data <- my.data[order(my.data$ID),]

n <- nrow(my.data) # Number of observations.
# Standardize the matrix of covariates.
Z <- scale(data.matrix(my.data[,c("Age", "PD.L1.MPS", "Hemoglobin")]))


# Set objects to save computing times.
ks.time <- rep(NA, max.sim)
caro.time <- rep(NA, max.sim)
nt.time <- rep(NA, max.sim)
ps.time <- rep(NA, max.sim)

# Ma and Hu (2013). Randomization algorithm based on Kernel Density.============
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
for (i in 1:max.sim){
  treat.assign <- kernel.randomization(i, Z = Z, n.zero = n.zero, bc.prob = bc.prob)
  ks.time[i] <- treat.assign$MeanTime
}


# CA-RO Algorithm of Bertsimas et al. (2019).===================================
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
for (i in 1:max.sim){
  treat.assign <- ca.ro(i, Z = Z, rho.p = rho.p,
                        upper.gamma.bound = upper.gamma.bound,
                        n.zero = n.zero, 
                        bc.prob = 1)
  caro.time[i] <- treat.assign$MeanTime
}


# Minimization algorithm for mean and variances of Nishi and Takaichi (2003).===
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
for (i in 1:max.sim){
  treat.assign <- minimization.mv(i, Z = Z, n.zero = n.zero,
                                  bc.prob = bc.prob)
  nt.time[i] <- treat.assign$MeanTime
}


# Pocock and Simon method.======================================================
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
for (i in 1:max.sim){
  treat.assign <- PS.minimization(i, Z = Z, n.categories = n.categories,
                                  n.zero = n.zero, 
                                  bc.prob = bc.prob)
  ps.time[i] <- treat.assign$MeanTime
}


# Save all designs.=============================================================
save(ks.time, caro.time, nt.time, ps.time,
     file = "PembroDesigns_CPUTime_Construction.RData")
