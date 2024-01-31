# ===================================================================
# Study the effect of rho in the CARO algorithm of Bertsimas et al. 
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
upper.gamma.bound <- 4


################################################################################
#########                        PERMBRO TRIAL.                 ################
################################################################################


# Develop simulation protocol.
p <- 3 # Number of covariates.
my.data <- read.table("data/Pembrolizumab.txt", header = T)
my.data <- my.data[order(my.data$ID),]

n <- nrow(my.data) # Number of observations.
# Standardize the matrix of covariates.
Z <- scale(data.matrix(my.data[,c("Age", "PD.L1.MPS", "Hemoglobin")]))


# CA-RO Algorithm of Bertsimas et al. (2019).===================================
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
results_caro_rhoOne <- apply(matrix(1:max.sim), MARGIN = 1, 
                             FUN = "ca.ro",
                             Z = Z, rho.p = 1,
                             upper.gamma.bound = upper.gamma.bound,
                             n.zero = n.zero, 
                             bc.prob = 1)

results_caro_rhoThree <- apply(matrix(1:max.sim), MARGIN = 1, 
                             FUN = "ca.ro",
                             Z = Z, rho.p = 3,
                             upper.gamma.bound = upper.gamma.bound,
                             n.zero = n.zero, 
                             bc.prob = 1)

results_caro_rhoSix <- apply(matrix(1:max.sim), MARGIN = 1, 
                      FUN = "ca.ro",
                      Z = Z, rho.p = 6,
                      upper.gamma.bound = upper.gamma.bound,
                      n.zero = n.zero, 
                      bc.prob = 1) 

results_caro_rhoTen <- apply(matrix(1:max.sim), MARGIN = 1, 
                             FUN = "ca.ro",
                             Z = Z, rho.p = 10,
                             upper.gamma.bound = upper.gamma.bound,
                             n.zero = n.zero, 
                             bc.prob = 1) 



# Save all designs.=============================================================
save(results_caro_rhoOne, results_caro_rhoThree, results_caro_rhoSix, results_caro_rhoTen, 
     file = "Study_Tuning_Parameters/PembroDesigns_CARO_Different_rho.RData")


################################################################################
#########                        SPASMS TRIAL.                 ################
################################################################################
# Develop simulation protocol.
p <- 2 # Number of covariates.
my.data <- read.table("data/Infant_Spasms.txt", header = T)
my.data <- my.data[order(my.data$ID),]

n <- nrow(my.data) # Number of observations.
# Standardize the matrix of covariates.
Z <- scale(data.matrix(my.data[,c("Age", "ISFreq")]))

# CA-RO Algorithm of Bertsimas et al. (2019).===================================
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
results_caro_rhoOne <- apply(matrix(1:max.sim), MARGIN = 1, 
                             FUN = "ca.ro",
                             Z = Z, rho.p = 1,
                             upper.gamma.bound = upper.gamma.bound,
                             n.zero = n.zero, 
                             bc.prob = 1)

results_caro_rhoThree <- apply(matrix(1:max.sim), MARGIN = 1, 
                               FUN = "ca.ro",
                               Z = Z, rho.p = 3,
                               upper.gamma.bound = upper.gamma.bound,
                               n.zero = n.zero, 
                               bc.prob = 1)

results_caro_rhoSix <- apply(matrix(1:max.sim), MARGIN = 1, 
                             FUN = "ca.ro",
                             Z = Z, rho.p = 6,
                             upper.gamma.bound = upper.gamma.bound,
                             n.zero = n.zero, 
                             bc.prob = 1) 

results_caro_rhoTen <- apply(matrix(1:max.sim), MARGIN = 1, 
                             FUN = "ca.ro",
                             Z = Z, rho.p = 10,
                             upper.gamma.bound = upper.gamma.bound,
                             n.zero = n.zero, 
                             bc.prob = 1) 


# Save all designs.=============================================================
save(results_caro_rhoOne, results_caro_rhoThree, results_caro_rhoSix, results_caro_rhoTen, 
     file = "Study_Tuning_Parameters/SpasmsDesigns_CARO_Different_rho.RData")

