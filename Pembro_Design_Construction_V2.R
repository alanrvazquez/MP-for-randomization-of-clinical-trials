# ===================================================================
# Comparisons of three adaptive randomization methods for
# continuous covariates and two treatments.
#
# Construction of designs using different strategies and parameter settings.
# 
# ===================================================================
library(randomizeR)
library(e1071) # For "moments" function.

source("code/adaptive_randomization.R")

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


# Ma and Hu (2013). Randomization algorithm based on Kernel Density.============
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
results_kr <- apply(matrix(1:max.sim), MARGIN = 1, 
                       FUN = "kernel.randomization",
                       Z = Z, n.zero = n.zero, 
                       bc.prob = bc.prob)



# CA-RO Algorithm of Bertsimas et al. (2019).===================================
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
results_caro <- apply(matrix(1:max.sim), MARGIN = 1, 
                         FUN = "ca.ro",
                         Z = Z, rho.p = rho.p,
                         upper.gamma.bound = upper.gamma.bound,
                         n.zero = n.zero, 
                         bc.prob = 1) 


# Random allocation.============================================================
set.seed(1389874)
results_random <- apply(matrix(1:max.sim), MARGIN = 1, 
                           FUN = "random", Z = Z, n.zero = n.zero) 


# Minimization algorithm for mean and variances of Nishi and Takaichi (2003).===
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
results_mv <- apply(matrix(1:max.sim), MARGIN = 1, 
                       FUN = "minimization.mv",
                       Z = Z, n.zero = n.zero,
                       bc.prob = bc.prob) 



# Pocock and Simon method.======================================================
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
results_ps<- apply(matrix(1:max.sim), MARGIN = 1, 
                      FUN = "PS.minimization",
                      Z = Z, n.categories = n.categories,
                      n.zero = n.zero, 
                      bc.prob = bc.prob) 


# Save all designs.=============================================================
save(results_kr, results_caro, results_mv, results_ps, results_random,
     file = "PembroDesigns_Individual_v2.RData")
